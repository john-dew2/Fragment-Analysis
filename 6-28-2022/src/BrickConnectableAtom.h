/*
 *  This file is part of esynth.
 *
 *  esynth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  esynth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with esynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _BRICK_CONNECTABLE_ATOM_GUARD
#define _BRICK_CONNECTABLE_ATOM_GUARD 1


#include <string>
#include <vector>


#include "AtomT.h"
#include "ConnectableAtom.h"


class Brick;


class BrickConnectableAtom : public ConnectableAtom
{
  protected:
    // Allowable atom types for connections
    AtomT** allowableTypes;

  public:
    bool CanConnectToAny() const { return false; }
    bool CanConnectTo(const Atom& that) const;

    BrickConnectableAtom(const BrickConnectableAtom* const that);
    BrickConnectableAtom(const std::string&, const Brick* const owner, const std::vector<std::string>& types);
    ~BrickConnectableAtom();

    bool IsLinkerAtom() const { return false; }
    bool IsBrickAtom() const { return true; }

    AtomT** getAllowableTypes() const { return allowableTypes; } 
    unsigned getNumAllowableTypes() const { return theAtom.numAllowConns; } 

    std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, BrickConnectableAtom& atom);
};

#endif
