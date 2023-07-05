/**
 *  OpenCCO implementation 
 *  Copyright (C) 2023 B. Kerautret;  Phuc Ngo, N. Passat H. Talbot and C. Jaquet
 *
 *  This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <https://www.gnu.org/licenses/>.
**/

#pragma once

#if defined(CONSTRUCTIONHELPERS_RECURSES)
#error Recursive header files inclusion detected in ConstructionHelpers.h
#else // defined(CONSTRUCTIONHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CONSTRUCTIONHELPERS_RECURSES

#if !defined CONSTRUCTIONHELPERS_h
/** Prevents repeated inclusion of headers. */
#define CONSTRUCTIONHELPERS_h

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/io/readers/GenericReader.h"
#include "CoronaryArteryTree.h"



namespace ExpandTreeHelpers {


template<typename DomCtr, int TDim>
inline
void
initFirstElemTree(CoronaryArteryTree< DomCtr, TDim > &aTree,
				  bool verbose = false)
{
	for(unsigned int i = 0; i < TDim; i++)
	{
		aTree.myTreeCenter[i] = static_cast<int>(aTree.myDomainController().myCenter[i]);
	}

	auto p = aTree.myDomainController().firstCandidatePoint();
	aTree.myVectSegments[0].myCoordinate = p;
}



template<typename DomCtr, int TDim>
inline
void
expandTree(CoronaryArteryTree< DomCtr, TDim > &aTree,
		   bool verbose = false, unsigned int nbMaxSearch = 100,
		   unsigned int nbTryCandidate = 100)
{
	srand ((unsigned int) time(NULL));
	bool isOK = false;
	unsigned int nbSeedFound = 1;
	unsigned int nbSeed = aTree.my_NTerm;

	for (unsigned int i = 1; i < nbSeed; i++)
	{
		DGtal::trace.progressBar(i, nbSeed);
		size_t nbSol = 0, itOpt = 0;
		CoronaryArteryTree< DomCtr, TDim > cTreeOpt = aTree;
		double volOpt = -1.0, vol = 0.0;
		unsigned int nbT = 0;

		while (nbSol==0 && nbT < nbMaxSearch)
		{
			auto pt = aTree.generateNewLocation(nbTryCandidate);
			nbT++;
			std::vector<unsigned int> vecN = aTree.getN_NearestSegments(pt,aTree.myNumNeighbor);
			for(size_t it=0; it<vecN.size(); it++)
			{
				auto ptBifurcation = aTree.findBarycenter(pt, vecN.at(it));
				if(!aTree.isIntersectingTree(pt, ptBifurcation,
											aTree.myVectSegments[vecN.at(it)].myRadius,
											vecN.at(it)))
				{
					CoronaryArteryTree< DomCtr, TDim  > cTree1 = aTree;
					isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor, verbose);
					if(isOK)
					{
						vol = cTree1.computeTotalVolume(1);
						
						if(volOpt<0.0 || volOpt > vol)
						{
							volOpt = vol;
							cTreeOpt = cTree1;
							itOpt = it;
						}

						nbSol++;
					}
				}
			}
		}

		if(nbT < nbMaxSearch)
		{
			nbSeedFound++;
		}

		aTree = cTreeOpt;
		aTree.updateLengthFactor();
		aTree.updateResistanceFromRoot();
		aTree.updateRootRadius();
	}

	if (nbSeed != nbSeedFound)
	{
		DGtal::trace.warning() << "All seeds not found due to too large domain constraints ("
							   << nbSeedFound << " over " << nbSeed << ")";
	}

	if (verbose)
	{
		std::cout << "====> Aperf=" << aTree.myRsupp*aTree.myRsupp*aTree.my_NTerm*M_PI
				  << " == " << aTree.my_aPerf << std::endl;
	}
}


/*
template<typename TTree>
inline
void
expandMultipleTrees(std::vector<TTree> & trees,
					bool verbose = false, unsigned int nbMaxSearch = 100,
					unsigned int nbTryCandidate = 100)
{
	srand ((unsigned int) time(NULL));
	
	std::vector<unsigned int> nbSeedFound(trees.size(), 1);
	std::vector<unsigned int> nbSeedAttempt(trees.size(), 1);

	std::vector<unsigned int> nbSeeds;
	for(auto it = trees.begin(); it != trees.end(); it++)
		nbSeeds.push_back(it->my_NTerm);

	for (unsigned int i = 1; i < nbSeed; i++)
	{
		DGtal::trace.progressBar(i, nbSeed);
		size_t nbSol = 0, itOpt = 0;
		TTree cTreeOpt = aTree;
		double volOpt = -1.0, vol = 0.0;
		unsigned int nbT = 0;

		while (nbSol==0 && nbT < nbMaxSearch)
		{
			auto pt = aTree.generateNewLocation(nbTryCandidate);
			nbT++;
			std::vector<unsigned int> vecN = aTree.getN_NearestSegments(pt,aTree.myNumNeighbor);
			for(size_t it=0; it<vecN.size(); it++)
			{
				auto ptBifurcation = aTree.findBarycenter(pt, vecN.at(it));
				if(!aTree.isIntersectingTree(pt, ptBifurcation,
											aTree.myVectSegments[vecN.at(it)].myRadius,
											vecN.at(it)))
				{
					TTree cTree1 = aTree;
					bool isOK = cTree1.isAddable(pt,vecN.at(it), 100, 0.01, cTree1.myNumNeighbor, verbose);
					if(isOK)
					{
						vol = cTree1.computeTotalVolume(1);
						// don't know why this code exists, removing it next build if everything works
						// if(volOpt<0.0)
						// {
						// 	volOpt = vol;
						// 	cTreeOpt = cTree1;
						// 	itOpt = it;
						// }
						// else
						// {
						// 	if(volOpt>vol)
						// 	{
						// 		volOpt = vol;
						// 		cTreeOpt = cTree1;
						// 		itOpt = it;
						// 	}
						// }

						volOpt = vol;
						cTreeOpt = cTree1;
						itOpt = it;
						nbSol++;
					}
				}
			}
		}

		if(nbT < nbMaxSearch)
		{
		nbSeedFound++;
		}

		aTree = cTreeOpt;
		aTree.updateLengthFactor();
		aTree.updateResistanceFromRoot();
		aTree.updateRootRadius();
	}

	if (nbSeed != nbSeedFound)
	{
		DGtal::trace.warning() << "All seeds not found due to too large domain constraints ("
		<< nbSeedFound << " over " << nbSeed << ")";
	}

	if (verbose){
		std::cout<<"====> Aperf="<<aTree.myRsupp*aTree.myRsupp*aTree.my_NTerm*M_PI<<" == "<<aTree.my_aPerf<<std::endl;
	}
}



template<typename TImage >
std::vector<std::vector<typename TImage::Domain::Point > >
getImageContours(const TImage &image,
unsigned int threshold=128){
std::vector<std::vector<typename TImage::Domain::Point > > v;
DGtal::trace.error() << "Use CCO is only implemented in 2D and 3D and use ImageContainerBySTLVector with unsigned char"
<< "to export domain."
<< "You use such an image : " << image <<  std::endl;
throw 1;
return v;
}


std::vector<std::vector<DGtal::Z2i::Point > >
getImageContours(const DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> &image,
unsigned int threshold){
typedef  DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> TImage;
typedef DGtal::functors::IntervalThresholder<typename TImage::Value> Binarizer;
DGtal::Z2i::KSpace ks;
if(! ks.init( image.domain().lowerBound(),
image.domain().upperBound(), true )){
DGtal::trace.error() << "Problem in KSpace initialisation"<< std::endl;
}

Binarizer b(threshold, 255);
DGtal::functors::PointFunctorPredicate<TImage,Binarizer> predicate(image, b);
DGtal::trace.info() << "DGtal contour extraction from thresholds ["<<  threshold << "," << 255 << "]" ;
DGtal::SurfelAdjacency<2> sAdj( true );
std::vector< std::vector< DGtal::Z2i::Point >  >  vectContoursBdryPointels;
DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAllPointContours4C( vectContoursBdryPointels, ks, predicate, sAdj );
return vectContoursBdryPointels;
}
*/



/**
Template Specialisation in 3D  to export surfel image  border of the restricted image domain.
* Todo @BK
*/
std::vector<std::vector<DGtal::Z3i::Point > >
getImageContours(const DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> &image,
unsigned int threshold){
std::vector<std::vector<DGtal::Z3i::Point > > v;

return v;
}



}

#endif // !defined CONSTRUCTIONHELPERS_h

#undef CONSTRUCTIONHELPERS_RECURSES
#endif // else defined(CONSTRUCTIONHELPERS_RECURSES)




