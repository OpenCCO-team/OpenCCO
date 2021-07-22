#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"




#pragma once

#if defined(CONSTRUCTIONHELPERS_RECURSES)
#error Recursive header files inclusion detected in ConstructionHelpers.h
#else // defined(CONSTRUCTIONHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CONSTRUCTIONHELPERS_RECURSES

#if !defined CONSTRUCTIONHELPERS_h
/** Prevents repeated inclusion of headers. */
#define CONSTRUCTIONHELPERS_h

namespace ConstructionHelpers {

void constructTree(double aPerf, int nbTerm,
                   std::string imageOrgan, bool verbose);

}

#endif // !defined CONSTRUCTIONHELPERS_h

#undef CONSTRUCTIONHELPERS_RECURSES
#endif // else defined(CONSTRUCTIONHELPERS_RECURSES)




