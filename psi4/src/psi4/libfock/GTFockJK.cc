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
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"

#include <sstream>
#include <vector>
#include <stdio.h>
#include <string>
#include <memory>


#include "psi4/libpsi4util/exception.h"
#include "jk.h"
#include "psi4/libgtfock_interface/Minimal_Interface.h"

using namespace psi;


namespace psi {
GTFock_JK::GTFock_JK(std::shared_ptr<BasisSet> Primary, size_t NMats, bool AreSymm) : JK(Primary)  {
//    common_init();
    /* Temporary containers. Data will be deep copied into GTFock */
    size_t count = 0UL;
    std::vector<int> Charges(primary_->molecule()->natom());
    std::vector<double> Places_X(primary_->molecule()->natom());
    std::vector<double> Places_Y(primary_->molecule()->natom());
    std::vector<double> Places_Z(primary_->molecule()->natom());
    std::vector<int> atom_shells(primary_->molecule()->natom());
    std::vector<int>  shell_prims(primary_->nshell());
    std::vector<int> ang_momentum(primary_->nshell());
    std::vector<double> coefficients(primary_->nprimitive());
    std::vector<double> exponents(primary_->nprimitive());

    /* int because natom is small */
    for (int i = 0; i < primary_->molecule()->natom(); i++ ) {
        Charges[i] = primary_->molecule()->charge(i);
        Places_X[i] = primary_->molecule()->x(i);
        Places_Y[i] = primary_->molecule()->y(i);
        Places_Z[i] = primary_->molecule()->z(i);
        atom_shells[i] = primary_->nshell_on_center(i);
    }
    for (size_t i = 0UL; i < primary_->nshell(); i++) {
        shell_prims[i] = primary_->shell(i).nprimitive();
        ang_momentum[i] = primary_->shell(i).am();
    // I would use memcpy, but I would have to dynamically cast to void*... 
    //    memcpy(&primary_->shell(i).coefs()[0], &coefficients[count], shell_prims[i]);
    //    memcpy(&primary_->shell(i).exps()[0], &exponents[count], shell_prims[i]);
        for(size_t j = 0UL; j < primary_->shell(i).nprimitive(); j++) {
            coefficients[count + j] = primary_->shell(i).coefs()[j];
            exponents[count + j] = primary_->shell(i).exps()[j];
        }
        count += shell_prims[i];
    }
//    for (size_t i = 0UL; i < primary_->nprimitive() ; i++ ){
//        coefficients[i] = primary_->
//        exponents[i] = 
//    }

    gtfock_ = std::make_shared<gtfock_interface>(NMats, AreSymm,
          primary_->molecule()->natom(), &Charges[0], 
          &Places_X[0], &Places_Y[0], &Places_Z[0],
          primary_->nprimitive(), primary_->nshell(), ( primary_->has_puream() ? 1 : 0 ),  
          &atom_shells[0], &shell_prims[0],
          &ang_momentum[0], &coefficients[0], &exponents[0] );
} /* GTFock_JK::GTFock_JK */
void GTFock_JK::compute_JK() {}
void GTFock_JK::print_header() const { outfile->Printf("Hi!\n"); }
GTFock_JK::~GTFock_JK() {}
} /* namespace psi */


#ifdef ENABLE_GTFOCK
#include <GTFock/MinimalInterface.h>
#else
namespace psi {
struct MinimalInterface {
    MinimalInterface(size_t, bool) { throw PSIEXCEPTION("PSI4 has not been compiled with GTFock support"); }
    void SetP(std::vector<SharedMatrix>&) {}
    void GetJ(std::vector<SharedMatrix>&) {}
    void GetK(std::vector<SharedMatrix>&) {}
};
}
#endif

#ifdef ENABLE_GTFOCK
namespace psi {
GTFockJK::GTFockJK(std::shared_ptr<psi::BasisSet> Primary) : JK(Primary), Impl_(new MinimalInterface()) {}
size_t GTFockJK::estimate_memory() {
    return 0; // Effectively
}
void GTFockJK::compute_JK() {
    NMats_ = C_left_.size();
    Impl_->create_pfock(NMats_, lr_symmetric_);
    Impl_->SetP(D_ao_);
    Impl_->GetJ(J_ao_);
    Impl_->GetK(K_ao_);
    Impl_->destroy_gtfock();
}
}
#endif
