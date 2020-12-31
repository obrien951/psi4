/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _GTFOCK_INTERFACE_H_
#define _GTFOCK_INTERFACE_H_

struct BasisSet;
struct PFock;
namespace psi {

class gtfock_interface {
    protected:
    BasisSet* gtf_basis;
    int NMat_; /* number of density matrices for fock calculations */
    bool AreSymm_;
    public:
    /* To deal with a naming conflict, the implementation does not see
 *     the psi basis. Instead, the data and metadata are supplied directly to
 *     the interface by the JK object that wraps it. 
 *     spa_{X,Y,Z} for "spatial {""} */
    gtfock_interface(int NMats, bool AreSymm,
                     int natoms, int *charges, 
                     double* spa_X, double* spa_Y, double* spa_Z,
                     int nprims, int nshells, int pure,
                     int *atom_shells, int *shell_prims,
                     int *ang_momentum, double *CC, double* alpha );
        
}; /* class gtfock_interface */
}  /* namespace psi */
#endif /*_GTFOCK_INTERFACE_H_ */
