#include "svg_elements.h"



namespace SVG
{

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               SVG::Animation                ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	// Constructors

	Animation::Animation()
	{
		//
	}


	Animation::Animation(const std::string & attribute_name,
					       const Timeline & timeline,
						   const bool repeat_count,		// = 0
						   const bool freeze)			// = false
				: m_attribute_name(attribute_name), m_repeat_count(repeat_count), m_freeze(freeze), m_timeline(timeline)
	{
		//
	}


	// Setters & Getters

	void Animation::setTimeline(const Timeline & timeline)
	{
		m_timeline = timeline;
	}


	Timeline & Animation::getTimelineRef()
	{
		return m_timeline;
	}


	std::string Animation::getAttributeName()
	{
		return m_attribute_name;
	}


	// Methods

	void Animation::print(std::ostream & os) const
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
////////////////////////             SVG::AnimatedElement            ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	// Constructors

	AnimatedElement::AnimatedElement()
		: m_animations()
	{
		//
	}

	// Methods

	void AnimatedElement::addAnimation(const Animation & animation)
	{
		m_animations.push_back(animation);
	}


	void AnimatedElement::printAnimations(std::ostream & os) const
	{
		for(const Animation & a : m_animations)
		{
			os << a;
		}
	}


	int AnimatedElement::animationCount() const
	{
		return m_animations.size();
	}


	Animation * AnimatedElement::getAnimation(const std::string & attribute_name)
	{
		for(Animation & animation : m_animations)
		{
			if(animation.getAttributeName() == attribute_name)		// found the corresponding animation
			{
				return &animation;
			}
		}

		return nullptr;			// didn't find it
	}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  SVG::Svg                   ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	// Constructors

	Svg::Svg()
		: m_x(), m_y(), m_width(), m_height()
	{
		//
	}


	Svg::Svg(double x, double y, double width, double height)
		: m_x(x), m_y(y), m_width(width), m_height(height)
	{

	}


	// Methods

	void Svg::print(std::ostream & os) const
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
		for(std::shared_ptr<AnimatedElement> ptr_e : m_contained_elements)
		{
			ptr_e->print(os);
		}

		// close the svg tag ( </svg> )
		os << "</svg>\n";
	}


	void Svg::addElement(const std::shared_ptr<AnimatedElement> & ptr_element)
	{
		m_contained_elements.push_back(ptr_element);
	}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  SVG::Rect                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	// Constructors

	Rect::Rect()
		: m_x(), m_y(), m_width(), m_height(), m_fill_color()
	{

	}


	Rect::Rect(double x, double y, double width, double height, const Color & fill_color)
		: m_x(x), m_y(y), m_width(width), m_height(height), m_fill_color(fill_color)
	{
		//
	}


	// Methods

	void Rect::print(std::ostream & os) const
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
			os << "</rect>\n";
		}
		else
		{
			os << "/>\n";
		}
	}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  SVG::Line                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



	// Constructors

	Line::Line()
		: m_x1(), m_y1(), m_x2(), m_y2(), m_thickness(), m_stroke_color()
	{
		//
	}


	Line::Line(double x1, double y1, double x2, double y2, double thickness, const Color & stroke_color)
		: m_x1(x1), m_y1(y1), m_x2(x2), m_y2(y2), m_thickness(thickness), m_stroke_color(stroke_color)
	{
		//
	}


	// Methods

	void Line::print(std::ostream & os) const
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
			os << "</line>\n";
		}
		else
		{
			os << "/>\n";
		}
	}
}