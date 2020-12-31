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

//extern "C" {
#include "Minimal_Interface.h"
#include "psi4/libgtfock_interface/GTFock/gtfock/pfock/pfock.h"
#include "psi4/libgtfock_interface/GTFock/libcint/CInt.h"
#include "psi4/libgtfock_interface/GTFock/GTMatrix/GTMatrix.h"
//}

psi::gtfock_interface::gtfock_interface(int NMats, bool AreSymm,
                                        int natoms, int *charges, 
                                        double* spa_X, double* spa_Y, double* spa_Z,
                                        int nprims, int nshells, int pure,
                                        int *atom_shells, int *shell_prims,
                                        int *ang_momentum, double *CC, double* alpha ){
    NMat_ = NMats;
    AreSymm_ = AreSymm;
    CInt_createBasisSet(&gtf_basis);
    CInt_importBasisSet(gtf_basis,
                        natoms, charges,
                        spa_X, spa_Y, spa_Z,
                        nprims, nshells, pure,
                        atom_shells,
                        shell_prims,
                        ang_momentum, CC, alpha);
}
