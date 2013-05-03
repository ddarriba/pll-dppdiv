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

#include "Alignment.h"
#include "MbRandom.h"
#include "Model.h"
#include "Parameter.h"
#include "Parameter_basefreq.h"
#include "Parameter_exchangeability.h"
#include "Parameter_expcalib.h"
#include "Parameter_rate.h"
#include "Parameter_shape.h"
#include "Parameter_speciaton.h"
#include "Parameter_tree.h"
#include "Parameter_cphyperp.h"
#include "Parameter_treescale.h"
#include "Calibration.h"
#include <string>
#include <vector>
#include <fstream>

using namespace std;

Model::Model(MbRandom *rp, Alignment *ap, string ts, double pm, double ra,
		double rb, double hal, double hbe, bool ubl, bool alnm, int offmv,
		bool rndNo, string clfn, int nodpr, double bdr, double bda,
		double fxclkrt, bool roofix, bool sfb, bool ehpc, bool dphpc,
		int dphpng, bool gamhp, int rmod, bool fxmod, tree *tr[2], DataType dt,
		int pmodel) {

	dataType = dt;
	proteinModel = pmodel;
	ranPtr = rp;
	alignmentPtr = ap;
	priorMeanN = pm;
	ranPtr->getSeed(startS1, startS2);
	runUnderPrior = false;
	treeTimePrior = nodpr;
	myCurLnL = 0.0;
	lnLGood = false;
	calibfilen = clfn;
	fixRootHeight = roofix;
	zeroNodeTimeMove = false;
	exponCalibHyperParm = ehpc;
	exponDPMCalibHyperParm = dphpc;
	turnedOffMove = offmv;
	fixedClockRate = fxclkrt;
	fixSomeModParams = fxmod;
	activeParm = 0;
	double initRootH = 1.0;
	rHtY = 0.0;
	rHtO = 0.0;
	int tsPrDist = 1;
	rootNExpRate = -1.0;
	bool rtCalib = false;
	pll_tree = tr;

	if (calibfilen.empty() == false) {
		initRootH = readCalibFile();
		Calibration *rCal = getRootCalibration();
		if (rCal != NULL && fixRootHeight == false) {
			rHtY = rCal->getYngTime();
			rHtO = rCal->getOldTime();
			tsPrDist = rCal->getPriorDistributionType();
			rtCalib = true;
			rootNExpRate = rCal->getCalExponRate();
		} else if (rCal == NULL && fixRootHeight == false) {
			rHtY = rHtY / 1.1;
			rHtO = -1.0;
			tsPrDist = 2;
		}
	}

	cout << "\nStarting with seeds: { " << startS1 << " , " << startS2
			<< " } \n\n";

	numGammaCats = 4;
	numPatterns = alignmentPtr->getNumChar();

	cpfix = false;
	if (turnedOffMove == 5)
		cpfix = true;
	else if (turnedOffMove == 6)
		cpfix = true;
	if (rmod > 1)
		cpfix = true;

	int nn = 2 * alignmentPtr->getNumTaxa() - 1;

	if (pm > nn - 1) {
		cerr << "ERROR: the prior on the mean number of tables (" << pm
				<< ") cannot exceed the number of nodes in the tree (" << nn
				<< ")!" << endl;
		exit(1);
	}

	Cphyperp *conp = new Cphyperp(ranPtr, this, hal, hbe, nn, pm, cpfix);
	ExpCalib *excal = new ExpCalib(ranPtr, this, dphpc, dphpng, initRootH,
			gamhp);
	NodeRate *nr = new NodeRate(ranPtr, this, nn, ra, rb, conp->getCurrentCP(),
			fxclkrt, rmod);
	for (int i = 0; i < 2; i++) {
		parms[i].push_back(new Basefreq(ranPtr, this, 4, fxmod, dataType, proteinModel)); // base frequency parameter
		parms[i].push_back(new Exchangeability(ranPtr, this)); // rate parameters of the GTR model
		parms[i].push_back(new Shape(ranPtr, this, numGammaCats, 2.0, fxmod)); // gamma shape parameter for rate variation across sites
		parms[i].push_back(
				new Tree(ranPtr, this, alignmentPtr, ts, ubl, alnm, rndNo,
						calibrs, initRootH, sfb, ehpc, excal, pll_tree[i])); // rooted phylogenetic tree
		parms[i].push_back(nr); // restaurant containing node rates
		parms[i].push_back(conp); // hyper prior on DPP concentration parameter
		parms[i].push_back(
				new Treescale(ranPtr, this, initRootH, rHtY, rHtO, tsPrDist,
						rtCalib, ehpc)); // the tree scale prior
		parms[i].push_back(new Speciation(ranPtr, this, bdr, bda)); // hyper prior on diversification for cBDP speciation
		parms[i].push_back(excal); // hyper prior exponential node calibration parameters
		//rearrangeModelParameters();
		switchActiveParm();
	}

	numParms = parms[0].size();

	for (int i = 0; i < numParms; i++)
		*parms[0][i] = *parms[1][i];

	for (int i = 0; i < numParms; i++)
		parms[0][i]->print(std::cout);

	setUpdateProbabilities(true);
	if (ehpc)
		excal->getAllExpHPCalibratedNodes();

	switchActiveParm();
	rearrangeModelParameters();
	switchActiveParm();
	rearrangeModelParameters();
	myCurLnL = lnLikelihood();
	cout << "lnL = " << myCurLnL << endl;
}

Model::~Model(void) {

}

void Model::switchActiveParm(void) {
#ifdef _USE_PTHREADS
	pll_tree[activeParm]->isActive = false;
#endif
	(activeParm == 0 ? activeParm = 1 : activeParm = 0);
#ifdef _USE_PTHREADS
	pll_tree[activeParm]->isActive = true;
#endif
}

void Model::switchActiveParm(int newActiveParm) {
	if (newActiveParm != activeParm)
		switchActiveParm();
}

Basefreq* Model::getActiveBasefreq(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Basefreq *derivedPtr = dynamic_cast<Basefreq *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

Tree* Model::getActiveTree(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Tree *derivedPtr = dynamic_cast<Tree *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

Treescale* Model::getActiveTreeScale(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Treescale *derivedPtr = dynamic_cast<Treescale *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

Exchangeability* Model::getActiveExchangeability(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Exchangeability *derivedPtr = dynamic_cast<Exchangeability *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

Shape* Model::getActiveShape(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Shape *derivedPtr = dynamic_cast<Shape *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

NodeRate* Model::getActiveNodeRate(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		NodeRate *derivedPtr = dynamic_cast<NodeRate *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

Speciation* Model::getActiveSpeciation(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Speciation *derivedPtr = dynamic_cast<Speciation *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

Cphyperp* Model::getActiveCphyperp(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		Cphyperp *derivedPtr = dynamic_cast<Cphyperp *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

ExpCalib* Model::getActiveExpCalib(void) {

	for (int i = 0; i < numParms; i++) {
		Parameter *p = parms[activeParm][i];
		ExpCalib *derivedPtr = dynamic_cast<ExpCalib *>(p);
		if (derivedPtr != 0)
			return derivedPtr;
	}
	return NULL;
}

void Model::checkModelParameters(void) {
	cout.precision(20);
	cout << "SHAPE" << endl;
	getActiveShape()->print(cout);
	if (dataType == NUCLEIC) {
		cout << "Frequencies" << endl;
		getActiveBasefreq()->print(cout);
		cout << "Rates" << endl;
		getActiveExchangeability()->print(cout);
	}
	double avg = 0.0;
	double max = 0.0;
	double min = 1.0;
	NodeRate *r = getActiveNodeRate();
	for (int i = 0; i < getActiveTree()->getNumNodes(); i++) {
		Node *p = getActiveTree()->getNodeByIndex(i);
		if (p->getAnc() != NULL) {
			cout << i << "\t";
			double v;
			if (p->getAnc()->getAnc() == NULL) {
				/* Root child */
				double branchProportionL = p->getAnc()->getNodeDepth()
						- p->getAnc()->getLft()->getNodeDepth();
				double rPL = r->getRateForNodeIndexed(
						p->getAnc()->getLft()->getIdx());
				double branchProportionR = p->getAnc()->getNodeDepth()
						- p->getAnc()->getRht()->getNodeDepth();
				double rPR = r->getRateForNodeIndexed(
						p->getAnc()->getRht()->getIdx());
				v = branchProportionL * rPL + branchProportionR * rPR;
			} else {
				double branchProportion = p->getAnc()->getNodeDepth()
						- p->getNodeDepth();
				double rP = r->getRateForNodeIndexed(p->getIdx());
				v = branchProportion * rP;
			}
			cout << v << endl;
			double bl_from_pll = -log(p->getBranchLength())
					* getActivePllTree()->fracchange;
			double diff = abs(bl_from_pll - v);
			if (diff > max) {
				max = diff;
			} else if (diff < min) {
				min = diff;
			}
			avg += diff;
			//cout.precision(20);
			//cout << "DIFF " << setw(4) << i << " " << bl_from_pll - v << endl;
		}
	}
	avg /= getActiveTree()->getNumNodes() - 1;
	cout << "Min: " << log10(min) << " Max: " << log10(max) << " Avg: "
			<< log10(avg) << endl;
}

double Model::lnLikelihood(bool fullTraversal) {
	if (runUnderPrior) {
		myCurLnL = 0.0;
	} else {
		if (fullTraversal)
			getActivePllTree()->start = getActivePllTree()->nodep[1];
		evaluateGeneric(getActivePllTree(), getActivePllTree()->start,
				fullTraversal);
		myCurLnL = getActivePllTree()->likelihood;
	}
	return myCurLnL;
}

Parameter* Model::pickParmToUpdate(int *id) {

	double u = ranPtr->uniformRv();
	double sum = 0.0;
	Parameter *parm = NULL;
	for (unsigned i = 0; i < updateProb.size(); i++) {
		sum += updateProb[i];
		if (u < sum) {
			*id = i;
			parm = parms[activeParm][i];
			break;
		}
	}
	return parm;
}

double Model::safeExponentiation(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}

void Model::rearrangeModelParameters(void) {
	initReversibleGTR(getActivePllTree(), CURRENT_PARTITION);
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
	masterBarrier(THREAD_COPY_INIT_MODEL, getActivePllTree());
#endif
	setPllBranchLengths();
}

void Model::setPllBranchLengths(void) {
	Tree *t = getActiveTree();
	NodeRate *r = getActiveNodeRate();
	for (int n = 0; n < t->getNumNodes(); n++) {
		Node *p = t->getDownPassNode(n);
		if (!IS_ROOT(p)) {
			setPllBranchLengths(p, r);
		}
	}
}

double Model::convertBranchLengthToPll(double branch_length) {
	return (exp(-branch_length / getActivePllTree()->fracchange));
}

void Model::setPllBranchLengths(Node *p, NodeRate *r) {
	int idx = p->getIdx();
	double prev = p->getBranchLength();
	double branchProportion = (p->getAnc()->getNodeDepth() - p->getNodeDepth()); // * getActiveTree()->getTreeScale();
	double rP = r->getRateForNodeIndexed(p->getIdx());
	double v = branchProportion * rP;
	double unrootedBranchLength;
	if (p->getAnc()->getAnc() == NULL) {
		Node *otherChild =
				p->getAnc()->getLft() == p ?
						p->getAnc()->getRht() : p->getAnc()->getLft();
		unrootedBranchLength = v
				+ (otherChild->getAnc()->getNodeDepth()
						- otherChild->getNodeDepth())
						* r->getRateForNodeIndexed(otherChild->getIdx());
	} else {
		unrootedBranchLength = v;
	}

	if (rP == 0.0) {
		cerr << "ERROR: Problem rP = 0" << endl;
		exit(1);
	}
	p->setBranchLength(convertBranchLengthToPll(unrootedBranchLength));
}

void Model::setTraversalDescriptor(Node *p) {
	setPllBranchLengths(p, getActiveNodeRate());
	if (p->getLft())
		setPllBranchLengths(p->getLft(), getActiveNodeRate());
	if (p->getRht())
		setPllBranchLengths(p->getRht(), getActiveNodeRate());
	if (!runUnderPrior) {
		newviewGeneric(getActivePllTree(), p->getPllNode(), false);
		newviewGeneric(getActivePllTree(), p->getPllNode()->back, false);
	}
	getActivePllTree()->start = p->getPllNode();
}

void Model::setSubtreePointers(Node *p) {
	if (!p->getIsLeaf()) {
		p->getPllNode()->x = 1;
		p->getPllNode()->next->x = 0;
		p->getPllNode()->next->next->x = 0;
		if (p->getLft()->getNodeDepth() <= p->getRht()->getNodeDepth()) {
			setSubtreePointers(p->getLft());
		} else {
			setSubtreePointers(p->getRht());
		}
	} else {
		getActivePllTree()->start = p->getPllNode();
	}
}

void Model::setNodeRateGrpIndxs(void) {

	Tree *t = getActiveTree();
	NodeRate *r = getActiveNodeRate();
	int rtID = t->getRoot()->getIdx();
	for (int n = 0; n < t->getNumNodes(); n++) {
		if (n != rtID) {
			Node *p = t->getNodeByIndex(n);
			int tn = r->getTableNumForNodeIndexed(n);
			double srt = r->getRateForNodeIndexed(n);
			p->setRtGrpIdx(tn);
			p->setRtGrpVal(srt);
		}
	}
}

void Model::updateAccepted(void) {

	int from = activeParm, to;
	if (from == 0)
		to = 1;
	else
		to = 0;
	for (int i = 0; i < numParms; i++)
		*parms[to][i] = *parms[from][i];

	switchActiveParm(to);
	rearrangeModelParameters();
	switchActiveParm(from);
}

void Model::updateRejected(void) {

	int to = activeParm, from;
	if (to == 0)
		from = 1;
	else
		from = 0;
	for (int i = 0; i < numParms; i++)
		*parms[to][i] = *parms[from][i];
	rearrangeModelParameters();
}

double Model::getMyCurrLnL(void) {

	if (lnLGood) {
		lnLGood = false;
		return myCurLnL;
	} else
		return lnLikelihood(true);

}

double Model::readCalibFile(void) {

	cout << "\nCalibrations:" << endl;
	bool rootIs = false;
	Calibration *rooCal;
	string ln = getLineFromFile(calibfilen, 1);
	int nlins = atoi(ln.c_str());
	int nnodes = alignmentPtr->getNumTaxa() - 1;
	string *calList = new string[nlins];
	for (int i = 0; i < nlins; i++) {
		calList[i] = getLineFromFile(calibfilen, i + 2);
		Calibration *cal = new Calibration(calList[i]);
		calibrs.push_back(cal);
		if (cal->getIsRootCalib()) {
			rooCal = cal;
			rootIs = true;
		}
	}
	delete[] calList;

	double initTScale = 1.0;
	double yb = 0.0;
	double ob = 0.0;
	if (rootIs) {
		if (rooCal->getPriorDistributionType() == 1) {
			yb = rooCal->getYngTime();
			ob = rooCal->getOldTime();
			if (yb == ob) {
				initTScale = yb;
				fixRootHeight = true;
			} else {
				initTScale = yb + (ranPtr->uniformRv() * (ob - yb));
				fixRootHeight = false;
			}
		} else if (rooCal->getPriorDistributionType() == 2) {
			fixRootHeight = false;
			yb = rooCal->getYngTime();
			double expMean = yb * 0.2;
			initTScale = yb + ranPtr->exponentialRv(1 / expMean);
		}
	} else {
		yb = 0.0;
		fixRootHeight = false;
		for (vector<Calibration *>::iterator v = calibrs.begin();
				v != calibrs.end(); v++) {
			double tmpv;
			if ((*v)->getPriorDistributionType() == 1)
				tmpv = (*v)->getOldTime();
			else if ((*v)->getPriorDistributionType() == 2)
				tmpv = (*v)->getYngTime() * 1.1;
			if (tmpv > yb)
				yb = tmpv;
		}
		ob = yb + (yb * 2);
		double tsc = yb + (ranPtr->uniformRv() * (ob - yb));
		initTScale = tsc;
		rHtY = yb;
		rHtO = ob;
	}

	if (nlins == nnodes) {
		bool fixall = true;
		for (vector<Calibration *>::iterator v = calibrs.begin();
				v != calibrs.end(); v++) {
			double tmpo = (*v)->getOldTime();
			double tmpy = (*v)->getYngTime();
			if (tmpo != tmpy) {
				fixall = false;
				break;
			}
		}
		zeroNodeTimeMove = fixall;
	} else if (nlins == nnodes - 1 && rootIs == false) {
		bool fixall = true;
		for (vector<Calibration *>::iterator v = calibrs.begin();
				v != calibrs.end(); v++) {
			double tmpo = (*v)->getOldTime();
			double tmpy = (*v)->getYngTime();
			if (tmpo != tmpy) {
				fixall = false;
				break;
			}
		}
		zeroNodeTimeMove = fixall;
	}

	cout << "\nInitial root height : " << initTScale << " [" << yb << ", " << ob
			<< "]" << endl;
	return initTScale;
}

Calibration* Model::getRootCalibration(void) {

	for (vector<Calibration *>::iterator v = calibrs.begin();
			v != calibrs.end(); v++) {
		if ((*v)->getIsRootCalib())
			return (*v);
	}
	return NULL;
}

void Model::setUpdateProbabilities(bool initial) {

	double bfp, srp, shp, ntp, dpp, cpa, tsp, spp, ehp;
	if (initial) {
		bfp = 0.3;
		srp = 0.3;
		shp = 0.3;
		ntp = 0.4;
		dpp = 0.5;
		cpa = 0.3;
		tsp = 0.3;
		spp = 0.4;
		ehp = 0.0;
	} else {
		bfp = 0.2;
		srp = 0.2;
		shp = 0.2;
		ntp = 0.4;
		dpp = 0.5;
		cpa = 0.3;
		tsp = 0.4;
		spp = 0.4;
		ehp = 0.0;
	}
	if (turnedOffMove == 1)
		bfp = 0.0;
	else if (turnedOffMove == 2)
		srp = 0.0;
	else if (turnedOffMove == 3)
		shp = 0.0;
	else if (turnedOffMove == 4) {
		ntp = 0.0;
		tsp = 0.0;
		spp = 0.0;
		ehp = 0.0;
	} else if (turnedOffMove == 5) {
		dpp = 0.0;
		cpa = 0.0;
		setNodeRateGrpIndxs();
		*parms[0][3] = *parms[1][3];
	} else if (turnedOffMove == 6 || cpfix == true)
		cpa = 0.0;
	else if (turnedOffMove == 8)
		spp = 0.0;
	if (fixRootHeight)
		tsp = 0.0;
	if (treeTimePrior == 1 || treeTimePrior == 4)
		spp = 0.0;
	if (zeroNodeTimeMove == 1) {
		ntp = 0.0;
		cout << "All internal node times are fixed" << endl;
	}
	if (fixedClockRate < 0.0)
		cpa = 0.0;
	if (exponCalibHyperParm) {
		if (initial)
			ehp = 0.3;
		else
			ehp = 0.4;
	}
	if (fixSomeModParams) {
		bfp = 0.0;
		shp = 0.0;
	}
	updateProb.clear();
	if (dataType == PROTEIC) {
		updateProb.push_back(0); // 1 basefreq
		updateProb.push_back(0); // 2 sub rates
	} else {
		updateProb.push_back(bfp); // 1 basefreq
		updateProb.push_back(srp); // 2 sub rates
	}
	updateProb.push_back(shp); // 3 gamma shape
	updateProb.push_back(ntp); // 4 node times
	updateProb.push_back(dpp); // 5 dpp rates
	updateProb.push_back(cpa); // 6 concentration parameter
	updateProb.push_back(tsp); // 7 tree scale parameter
	updateProb.push_back(spp); // 8 speciation parameters
	updateProb.push_back(ehp); // 9 exponential calibration hyper priors
	double sum = 0.0;
	for (unsigned i = 0; i < updateProb.size(); i++)
		sum += updateProb[i];
	for (unsigned i = 0; i < updateProb.size(); i++)
		updateProb[i] /= sum;
}

void Model::printDescription(ostream &o) {
	if (dataType == NUCLEIC) {
		getActiveBasefreq()->print(o);
		getActiveExchangeability()->print(o);
	}
	getActiveShape()->print(o);
	getActiveSpeciation()->print(o);
	getActiveNodeRate()->print(o);
	getActiveTreeScale()->print(o);
	getActiveTree()->print(o);
}
