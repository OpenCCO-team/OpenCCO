#include "CoronaryArteryTree.h"
#include "GeomHelpers.h"
#include "ExpandTreeHelpers.h"
#include "XmlHelpers.h"


int main()
{
	typedef SquareDomainCtrl<2> Dom;
	typedef CoronaryArteryTree<Dom, 2> TTree;
	Dom::TPoint pCenter(0,0);
	Dom d(1.0, pCenter);

	TTree tree(2000, 4, d, 1.0);

	ExpandTreeHelpers::initSecondRoot(tree);

	//ExpandTreeHelpers::expandTree(tree);

	//XMLHelpers::writeTreeToXml(tree, "tree_init_2D.xml");

	return 0;
}