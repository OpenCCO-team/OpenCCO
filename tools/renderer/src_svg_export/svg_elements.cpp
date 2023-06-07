#include "svg_elements.h"


/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 SVGAnimate                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



// Constructors

SVGAnimate::SVGAnimate()
{
	//
}


SVGAnimate::SVGAnimate(const std::string & attribute_name,
				       const Timeline & timeline,
					   const bool repeat_count,		// = 0
					   const bool freeze)			// = false
			: m_attribute_name(attribute_name), m_repeat_count(repeat_count), m_freeze(freeze), m_timeline(timeline)
{
	//
}


// Setters

void SVGAnimate::setTimeline(const Timeline & timeline)
{
	m_timeline = timeline;
}


// Methods

void SVGAnimate::print(std::ostream & os) const
{
	// tag
	os << "<animate ";
	
	// attributes
	streamAttribute(os, "attributeName", m_attribute_name);
	
	if(m_repeat_count <= 0) 
		streamAttribute(os, "repeatCount", "indefinite");
	else
		streamAttribute(os, "repeatCount", m_repeat_count);

	if(m_freeze) { streamAttribute(os, "fill", "freeze"); }

	os << m_timeline << "/>\n";
}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////             SVGAnimatedElement              ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



// Constructors

SVGAnimatedElement::SVGAnimatedElement()
	: m_animations()
{
	//
}

// Methods

void SVGAnimatedElement::addAnimation(const SVGAnimate & animation)
{
	m_animations.push_back(animation);
}


void SVGAnimatedElement::printAnimations(std::ostream & os) const
{
	for(const SVGAnimate & a : m_animations)
	{
		os << a;
	}
}


int SVGAnimatedElement::animationCount() const
{
	return m_animations.size();
}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   SVGSvg                    ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


// Constructors

SVGSvg::SVGSvg()
	: m_x(), m_y(), m_width(), m_height()
{
	//
}


SVGSvg::SVGSvg(double x, double y, double width, double height)
	: m_x(x), m_y(y), m_width(width), m_height(height)
{

}


// Methods

void SVGSvg::print(std::ostream & os) const
{
	// open the svg tag ( <svg ...> )
 	os  << "<svg ";
	streamAttribute(os, "width", 15, "cm");
	streamAttribute(os, "height", 15, "cm");

	os << "viewBox=\"" << m_x << " " << m_y << " " << m_width << " " << m_height << "\" ";

	streamAttribute(os, "xmlns", "http://www.w3.org/2000/svg");
	os << ">\n";

	printAnimations(os);

	// the SVG elements contained withing the svg tag ( <svg> ... </svg> )
	for(std::shared_ptr<SVGAnimatedElement> ptr_e : m_contained_elements)
	{
		ptr_e->print(os);
	}

	// close the svg tag ( </svg> )
	os << "</svg>\n";
}


void SVGSvg::addElement(const std::shared_ptr<SVGAnimatedElement> & ptr_element)
{
	m_contained_elements.push_back(ptr_element);
}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   SVGRect                   ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



// Constructors

SVGRect::SVGRect()
	: m_x(), m_y(), m_width(), m_height(), m_fill_color()
{

}


SVGRect::SVGRect(double x, double y, double width, double height, const Color & fill_color)
	: m_x(x), m_y(y), m_width(width), m_height(height), m_fill_color(fill_color)
{
	//
}


// Methods

void SVGRect::print(std::ostream & os) const
{
	os  << "<rect ";
	streamAttribute(os, "x", m_x);	
	streamAttribute(os, "y", m_y);
	streamAttribute(os, "width", m_width);
	streamAttribute(os, "height", m_height);
	streamAttribute(os, "fill", m_fill_color);

	if(animationCount() > 0)
	{
		os << ">\n";
		printAnimations(os);
		os << "</rect>";
	}
	else
	{
		os << "/>\n";
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   SVGLine                   ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



// Constructors

SVGLine::SVGLine()
	: m_x1(), m_y1(), m_x2(), m_y2(), m_thickness(), m_stroke_color()
{
	//
}


SVGLine::SVGLine(double x1, double y1, double x2, double y2, double thickness, const Color & stroke_color)
	: m_x1(x1), m_y1(y1), m_x2(x2), m_y2(y2), m_thickness(thickness), m_stroke_color(stroke_color)
{
	//
}


// Methods

void SVGLine::print(std::ostream & os) const
{
	os  << "<line ";

	streamAttribute(os, "x1", m_x1);	
	streamAttribute(os, "y1", m_y1);
	streamAttribute(os, "x2", m_x2);	
	streamAttribute(os, "y2", m_y2);
	streamAttribute(os, "stroke-width", m_thickness);
	streamAttribute(os, "stroke", m_stroke_color);
	streamAttribute(os, "style", "stroke-linecap:round;stroke-linejoin:miter");

	if(animationCount() > 0)
	{
		os << ">\n";
		printAnimations(os);
		os << "</line>";
	}
	else
	{
		os << "/>\n";
	}
}