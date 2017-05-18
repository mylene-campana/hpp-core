//
// Copyright (c) 2014 CNRS
// Authors: Florent Lamiraux
//
// This file is part of hpp-core
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#include <hpp/core/path-planner.hh>

#include <hpp/util/debug.hh>

#include <hpp/core/roadmap.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/problem-target.hh>
#include <hpp/core/node.hh>
#include <hpp/core/edge.hh>
#include <hpp/core/path.hh>
#include <hpp/core/path-validation.hh>
#include <hpp/core/path-projector.hh>
#include <hpp/core/steering-method.hh>
#include "astar.hh"

# include <hpp/core/steering-method-straight.hh>
#include <iomanip> // for std::setprecision

namespace hpp {
  namespace core {

    PathPlanner::PathPlanner (const Problem& problem) :
      problem_ (problem), roadmap_ (Roadmap::create (problem.distance (),
						     problem.robot())),
      interrupt_ (false)
    {
    }

    PathPlanner::PathPlanner (const Problem& problem,
			      const RoadmapPtr_t& roadmap) :
      problem_ (problem), roadmap_ (roadmap),
      interrupt_ (false)
    {
    }

    void PathPlanner::init (const PathPlannerWkPtr_t& weak)
    {
      weakPtr_ = weak;
    }

    const RoadmapPtr_t& PathPlanner::roadmap () const
    {
      return roadmap_;
    }

    const Problem& PathPlanner::problem () const
    {
      return problem_;
    }

    void PathPlanner::startSolve ()
    {
      problem_.checkProblem ();
      // Tag init and goal configurations in the roadmap
      roadmap()->initNode (problem_.initConfig ());
      problem_.target()->initRoadmap ();
    }

    PathVectorPtr_t PathPlanner::solve ()
    {
      /*value_type timmy = time (NULL);
	srand (timmy);
	hppDout (info, "time (NULL)= " << std::setprecision (15) << timmy);*/
      interrupt_ = false;
      std::size_t oneStepCount = 0;
      const std::size_t plannerIterLimit = problem_.plannerIterLimit ();
      hppDout (info, "plannerIterLimit= " << plannerIterLimit);
      bool solved = false;
      startSolve ();
      hppDout (info, "check if already solved");
      solved = roadmap()->pathExists ();
      if (solved ) {
	hppDout (info, "configs already connected");
      } else {
	hppDout (info, "configs not already connected, try direct path");
	tryDirectPath ();
	solved = roadmap()->pathExists ();
	if (solved ) {
	  hppDout (info, "tryDirectPath succeeded");
	}
      }
      if (interrupt_) throw std::runtime_error ("Interruption");
      while (!solved && oneStepCount < plannerIterLimit) {
        problem_.target ()->oneStep ();
	oneStep ();
	oneStepCount++;
	hppDout (info, "oneStepCount= " << oneStepCount);
	solved = roadmap()->pathExists ();
	if (interrupt_) throw std::runtime_error ("Interruption");
      }
      PathVectorPtr_t planned = PathVector::create (problem_.robot ()->configSize (), problem_.robot ()->numberDof ());
      if (oneStepCount < plannerIterLimit)
	planned = computePath ();
      else { // return directPath (even if invalid)
	problem_.nbPathPlannerFails_++;
	hppDout (info, "add (maybe invalid) directPath");
	const SteeringMethodPtr_t& sm (problem ().steeringMethod ());
	PathProjectorPtr_t pathProjector (problem ().pathProjector ());
	PathPtr_t projPath, path;
	NodePtr_t initNode = roadmap ()->initNode();
	for (Nodes_t::const_iterator itn = roadmap ()->goalNodes ().begin();
	     itn != roadmap ()->goalNodes ().end (); ++itn) {
	  ConfigurationPtr_t q1 ((initNode)->configuration ());
	  ConfigurationPtr_t q2 ((*itn)->configuration ());
	  hppDout (info, "before sm");
	  if (!sm)
	    throw std::runtime_error ("no steering-method in problem");
	  path = (*sm) (*q1, *q2);
	  hppDout (info, "before path projector");
	  if (pathProjector) {
	    hppDout (info, "try to project");
	    if (!pathProjector->apply (path, projPath)) continue;
	  } else {
	    hppDout (info, "try NOT to project");
	    if (path) {
	      hppDout (info, "found steering-method path");
	      projPath = path;
	    }
	    else {
	      hppDout (info, "empty steering-method path");
	      throw std::runtime_error ("empty steering-method path");
	    }
	  }
	  if (projPath) {
	    hppDout (info, "add edges");
            roadmap ()->addEdge (initNode, *itn, projPath);
	    roadmap ()->addEdge (*itn, initNode, projPath->reverse ());
	    planned->appendPath (projPath);
	  } else { // constraints not applied, but still return directPath
	    roadmap ()->addEdge (initNode, *itn, path);
	    roadmap ()->addEdge (*itn, initNode, path->reverse ());
	  }
	  planned->appendPath (path);
	}
      }
      hppDout (info, "before finishSolve");
      return finishSolve (planned);
    }

    void PathPlanner::interrupt ()
    {
      interrupt_ = true;
    }

    PathVectorPtr_t PathPlanner::computePath () const
    {
      Astar astar (roadmap(), problem_.distance ());
      return astar.solution ();
    }

    PathVectorPtr_t PathPlanner::finishSolve (const PathVectorPtr_t& path)
    {
      return path;
    }

    void PathPlanner::tryDirectPath ()
    {
      // call steering method here to build a direct conexion
      const SteeringMethodPtr_t& sm (problem ().steeringMethod ());
      PathValidationPtr_t pathValidation (problem ().pathValidation ());
      PathProjectorPtr_t pathProjector (problem ().pathProjector ());
      PathPtr_t validPath, projPath, path;
      NodePtr_t initNode = roadmap ()->initNode();
      for (Nodes_t::const_iterator itn = roadmap ()->goalNodes ().begin();
	   itn != roadmap ()->goalNodes ().end (); ++itn) {
	ConfigurationPtr_t q1 ((initNode)->configuration ());
	ConfigurationPtr_t q2 ((*itn)->configuration ());
	assert (*q1 != *q2);
	path = (*sm) (*q1, *q2);
        if (!path) continue;
        if (pathProjector) {
          if (!pathProjector->apply (path, projPath)) continue;
        } else {
          projPath = path;
        }
        if (projPath) {
	  PathValidationReportPtr_t report;
          bool pathValid = pathValidation->validate (projPath, false, validPath,
						     report);
          if (pathValid && validPath->timeRange ().second !=
              path->timeRange ().first) {
	    hppDout (info, "add corresponding edges in RM");
            roadmap ()->addEdge (initNode, *itn, projPath);
	    roadmap ()->addEdge (*itn, initNode, projPath->reverse ());
          }
        }
      }
    }
  } //   namespace core
} // namespace hpp
