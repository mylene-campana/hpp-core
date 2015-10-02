//
// Copyright (c) 2015 CNRS
// Authors: Florent Lamiraux, Mylene Campana
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

#ifndef HPP_CORE_PRUNE_HH
# define HPP_CORE_PRUNE_HH

# include <hpp/core/path-optimizer.hh>

namespace hpp {
  namespace core {
    /// \addtogroup path_optimization
    /// \{

    /// Prune
    ///
    /// Remove redundant nodes as proposed in:
    /// http://www.staff.science.uu.nl/~gerae101/motion_planning/clearance.html
    ///
    /// \note The optimizer assumes that the input path is a vector of optimal
    ///       paths for the distance function.
    class HPP_CORE_DLLAPI Prune : public PathOptimizer
    {
    public:
      /// Return shared pointer to new object.
      static PrunePtr_t create (const Problem& problem);

      /// Optimize path
      virtual PathVectorPtr_t optimize (const PathVectorPtr_t& path);
    protected:
      Prune (const Problem& problem);
    }; // class Prune
    /// \}
  } // namespace core
} // namespace hpp
#endif // HPP_CORE_PRUNE_HH
