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

#include "Parameter.h"
#include "Parameter_shape.h"
#include "Parameter_tree.h"
#include "MbRandom.h"
#include "Model.h"
#include "util.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

Shape::Shape(MbRandom *rp, Model *mp, int nc, double lam, bool fx) :
		Parameter(rp, mp, mp->getActivePllTree()) {

	numCats = nc;
	lambda = lam;
	alpha = &pll_tree->partitionData[CURRENT_PARTITION].alpha;
	double minA = ALPHA_MIN;
	double maxA = 300.0;
        double newAlpha;
	do {
		newAlpha = ranPtr->exponentialRv(lambda);
	} while (newAlpha < minA | newAlpha > maxA);
        *alpha = newAlpha;
	rates =	pll_tree->partitionData[CURRENT_PARTITION].gammaRates;
	makeGammaCats(*alpha, rates, numCats,
			pll_tree->useMedian);
	name = "SH";
	if (fx) {
		*alpha = 1.0;
	}
}

Shape::~Shape(void) {

}

Shape& Shape::operator=(const Shape &s) {

	if (this != &s)
		clone(s);
	return *this;
}

void Shape::clone(const Shape &s) {

	if (numCats == s.numCats) {
		lambda = s.lambda;
		*alpha = *s.alpha;
		for (int i = 0; i < numCats; i++)
			rates[i] = s.rates[i];
	} else {
		cerr << "ERROR: Expected gamma rate vectors to be of equal size."
				<< endl;
		exit(1);
	}
}

void Shape::print(std::ostream & o) const {

	o << "Gamma Shape: " << *alpha << "( ";
	for (int i = 0; i < numCats; i++)
		o << fixed << rates[i] << " ";
	o << ")" << endl;
}

double Shape::update(double &oldLnL) {
	double tuning = log(2.0);
	double oldAlpha = *alpha;
	double rv = ranPtr->uniformRv();
	double newAlpha = oldAlpha * exp(tuning * (rv - 0.5));

	bool validAlph = false;
	double minA = ALPHA_MIN;
	double maxA = 300.0;
	do {
		if (newAlpha < minA)
			newAlpha = minA * minA / newAlpha;
		else if (newAlpha > maxA)
			newAlpha = maxA * maxA / newAlpha;
		else
			validAlph = true;
	} while (!validAlph);

	updateGammaRateCats(newAlpha);
	*alpha = newAlpha;

	double lnProposalRatio = log(newAlpha) - log(oldAlpha);
	double lnPriorRatio = lambda * (oldAlpha - newAlpha);

	modelPtr->rearrangeModelParameters();
	return lnPriorRatio + lnProposalRatio;
}

double Shape::lnPrior(void) {

	return ranPtr->lnExponentialPdf(lambda, *alpha);
}

void Shape::updateGammaRateCats(double alph) {

	if (numCats == 1)
		rates[0] = 1.0;
	else
		makeGammaCats(alph, rates, numCats,
				pll_tree->useMedian);
}

string Shape::writeParam(void) {

	stringstream ss;
	ss << "Gamma Shape: " << *alpha << "( ";
	for (int i = 0; i < numCats; i++)
		ss << fixed << setprecision(4) << rates[i] << " ";
	ss << ")" << endl;
	string outp = ss.str();
	return outp;
}
