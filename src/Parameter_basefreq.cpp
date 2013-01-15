/* 
 * PLL-DPPDiv 1.0 source code
 * Copyright (c) 2012
 * Diego Darriba (1,2) (ddarriba@udc.es)
 * Andre J. Aberer (3) (Andre.Aberer@h-its.org)
 * Tomas Flouri (3) (Tomas.Flouri@h-its.org)
 * Fernando Izquierdo-Carrasco (3) (Fernando.Izquierdo@h-its.org)
 * Tracy Heath (4) (tracy@berkeley.edu)
 * Alexandros Stamatakis (3,5) (Alexandros.Stamatakis@h-its.org)
 *
 * (1) Electronics and Systems, University of A Coruña, Spain
 * (2) Biochemistry, Genetics and Inmunology, University of Vigo, Spain
 * (3) Heidelberg Institute for Theoretical Studies, Heidelberg, Germany
 * (4) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (5) Karlsruhe Institute of Technology. Institute for Theoretical Informatics, Karlsruhe, Germany.
 *
 * DPPDiv version 1.0b source code (git: 9c0ac3d2258f89827cfe9ba2b5038f0f656b82c1)
 * Copyright 2009-2011
 * Tracy Heath(1,2,3) (NSF postdoctoral fellowship in biological informatics DBI-0805631)
 * Mark Holder(1)
 * John Huelsenbeck(2)
 *
 * (1) Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS 66045
 * (2) Integrative Biology, University of California, Berkeley, CA 94720-3140
 * (3) email: tracyh@berkeley.edu
 *
 * DPPDiv is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * Some of this code is from publicly available source by John Huelsenbeck
 */

#include "MbRandom.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_tree.h"
#include "util.h"

#include <iostream>
#include <iomanip>

using namespace std;

Basefreq::Basefreq(MbRandom *rp, Model *mp, int ns, bool fx) : Parameter(rp, mp, mp->getActivePllTree()) {

	numStates = ns;
	freqs = pll_tree->partitionData[CURRENT_PARTITION].frequencies;
	alpha = MbVector<double>(numStates);
	for (int i=0; i<numStates; i++)
		alpha[i] = 1.0;
	ranPtr->dirichletRv(alpha, freqs);
	alpha0 = 500.0;
	name = "BF";
	if(fx){
		for (int i=0; i<numStates; i++)
			freqs[i] = 0.25;
	}
}

Basefreq::~Basefreq(void) {

}

Basefreq& Basefreq::operator=(const Basefreq &b) {

	if (this != &b)
		clone(b);
	return *this;
}

void Basefreq::clone(const Basefreq &b) {

	if (b.numStates == numStates)
		{
		for (int i=0; i<numStates; i++)
			freqs[i] = b.freqs[i];
		}
	else
		{
		cerr << "ERROR: Expected base frequency vectors of equal size." << endl;
		exit(1);
		}
}

void Basefreq::print(std::ostream & o) const {

	o << "Base Frequency: ";
	for (int i=0; i<numStates; i++)
		o << fixed << freqs[i] << " ";
	o << endl;
}

double Basefreq::update(double &oldLnL) {

	MbVector<double> aForward(numStates);
	MbVector<double> aReverse(numStates);
	MbVector<double> oldFreqs(numStates);
	MbVector<double> newFreqs(numStates);
	
	for (int i=0; i<numStates; i++)
		{
		oldFreqs[i] = freqs[i];
		aForward[i] = freqs[i] * alpha0;
		}
		
	ranPtr->dirichletRv(aForward, newFreqs);
	
	double sum = 0.0;
	for(int i=0; i<numStates; i++){
		if(newFreqs[i] < 0.0001)
			newFreqs[i] = 0.0001;
		sum += newFreqs[i];
	}
	for(int i=0; i<numStates; i++)
		newFreqs[i] /= sum;
	
	for (int i=0; i<numStates; i++)
		freqs[i] = newFreqs[i];
	
	for (int i=0; i<numStates; i++)
		aReverse[i] = newFreqs[i] * alpha0;
	
	double lnProposalRatio = ranPtr->lnDirichletPdf(aReverse, oldFreqs) - ranPtr->lnDirichletPdf(aForward, newFreqs);

	modelPtr->rearrangeModelParameters();

	return lnProposalRatio;
}

double Basefreq::lnPrior(void) {

	return ranPtr->lnDirichletPdf(alpha, freqs);
}

string Basefreq::writeParam(void){
	
	stringstream ss;
	ss << "Base Frequency: ";
	for (int i=0; i<numStates; i++)
		ss << fixed << setprecision(4) << freqs[i] << " ";
	ss << endl;
	string outp = ss.str();
	return outp;
}
