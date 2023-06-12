#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include "svg_attributes.h"



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 SVGAnimate                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief SVG animate element, inherits from Streamable
 * @brief Convoluted class structure, but the logic is that an <animate> tag 
 * @brief cannot have <animate> tags nested within
 **/
class SVGAnimate : public Streamable
{
public:
	// Constructors
	
	/**
	 * @brief Constructor, initializes all members with default constructor
	 **/
	SVGAnimate();

	/**
	 * @brief Constructor, initializes all members with parameters
	 **/
	SVGAnimate(const std::string & attribute_name,
		       const Timeline & timeline,
			   const bool repeat_count = 0,
			   const bool freeze = false);

	// Setters & Getters

	/**
	 * @brief Sets the timeline of the animation 
	 * @param timeline The timeline
	 **/
	void setTimeline(const Timeline & timeline);

	/**
	 * @brief getter with the purpose of editing the member variable m_timeline
	 * @returns a reference to m_timeline
	 **/
	Timeline & getTimelineRef();

	/**
	 * @brief getter to access the value of m_attribute_name
	 * @returns a reference to m_timeline
	 **/
	std::string getAttributeName();

	// Methods
	
	/**
	 * @brief Implementation of the pure virtual function from Streamable
	 * @brief Outputs the SVGAnimate to the stream parameter 
	 * @param os The output stream
	 **/
	void print(std::ostream & os) const override;

private:
	std::string m_attribute_name;		// Name of the animated attribute
	Timeline m_timeline;				// Animation duration and key times
	int m_repeat_count;					// Number of animation loops (value <= 0 will be treated as "indefinite")
	bool m_freeze;						// Whether the animation freezes at the end or snaps back to the initial state
};	



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////             SVGAnimatedElement              ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief SVG animated element, inherits from Streamable
 * @brief Abstract class, any SVG element that support animation will inherit from SVGAnimatedElement
 * @brief i.e elements that can have nested <animate> tags
 **/
class SVGAnimatedElement : public Streamable
{
public:
	// Constructors

	/**
	 * @brief Constructor, initializes all members with default constructor
	 **/
	SVGAnimatedElement();

	// Methods

	/**
	 * @brief Pure virtual method, output the SVGAnimatedElement to the stream 
	 * @brief (called by operator<<(std::ostream & os, const Streamable & element) )
	 * @param os The output stream
	 **/
	virtual void print(std::ostream & os) const = 0;

	/**
	 * @brief Appends an SVG animation to m_animations
	 * @brief Checks if an animation for the same attribute exists or not, and by default will overwrite it
	 * @param animation The animation to append
	 **/
	void addAnimation(const SVGAnimate & animation);

	/**
	 * @brief Outputs m_animations to the stream
	 * @param os The stream output to
	 **/
	void printAnimations(std::ostream & os) const;

	/**
	 * @brief Number of animations in m_animations
	 * @returns m_animation.size()
	 **/
	int animationCount() const;

	/**
	 * @brief Finds the animation associated with the specified attribute_name
	 * @returns a pointer to the SVGAnimate element if it exists, nullptr if it doesn't
	 **/
	SVGAnimate * getAnimation(const std::string & attribute_name);

private:
	std::vector<SVGAnimate> m_animations;
};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   SVGSvg                    ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief SVG svg, inherits from SVGAnimatedElement
 * @brief Represents the <svg> ... </svg> tag
 * @brief Also acts as a container of other svg elements
 **/
class SVGSvg : public SVGAnimatedElement
{
public:
	// Constructors

	/**
	 * @brief Constructor, initializes all members with default constructor
	 **/
	SVGSvg();

	/**
	 * @brief Constructor, initializes all members with parameters
	 * @param x The x coordinate of the top left corner of the SVG
	 * @param y The y coordinate of the top left corner of the SVG
	 * @param width The width of the SVG
	 * @param height The height of the SVG
	 **/
	SVGSvg(double x, double y, double width, double height);

	// Methods

	/**
	 * @brief Implementation of the pure virtual function from SVGAnimatedElement
	 * @brief Outputs the SVGSvg to the stream parameter 
	 * @param os The output stream
	 **/
	void print(std::ostream & os) const override;

	/**
	 * @brief Appends a shared pointer of an SVGAnimatedElement to m_contained_elements
	 * @param ptr_element The shared pointer to append
	 **/
	void addElement(const std::shared_ptr<SVGAnimatedElement> & ptr_element);

private:
	// Members
	double m_x;					// Top left corner x coordinate
	double m_y;					// Top left corner y coordinate
	double m_width;				// SVG width
	double m_height;			// SVG height
	std::vector< std::shared_ptr<SVGAnimatedElement> > m_contained_elements;	// svg elements within the svg tags
};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   SVGRect                   ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief SVG rectangle, inherits from SVGAnimatedElement
 * @brief Represents the <rect> ... </rect> tag
 **/
class SVGRect : public SVGAnimatedElement
{
public:
	// Constructors

	/**
	 * @brief Constructor, initializes all members with default constructor
	 **/
	SVGRect();

	/**
	 * @brief Constructor, initializes all members with parameters
	 * @param x The x coordinate of the top left corner of the rectangle
	 * @param y The y coordinate of the top left corner of the rectangle
	 * @param width The width of the rectangle
	 * @param height The height of the rectangle
	 * @param fill_color The color of the rectangle
	 **/
	SVGRect(double x, double y, double width, double height, const Color & fill_color);

	// Methods

	/**
	 * @brief Implementation of the pure virtual function from SVGAnimatedElement
	 * @brief Outputs the SVGRect to the stream parameter 
	 * @param os The output stream
	 **/
	void print(std::ostream & os) const override;

private:
	// Members
	double m_x;					// Top left corner x coordinate
	double m_y;					// Top left corner y coordinate
	double m_width;				// Rectangle width
	double m_height;			// Rectangle height
	Color m_fill_color;			// RGB color of the rectangle
};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   SVGLine                   ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief SVG line, inherits from SVGAnimatedElement
 * @brief Represents the <line> ... </line> tag
 **/
class SVGLine : public SVGAnimatedElement
{
public:
	// Constructors

	/**
	 * @brief Constructor, initializes all members with default constructor
	 **/
	SVGLine();

	/**
	 * @brief Constructor, initializes all members with parameters
	 * @param x1 The x coordinate of the first point
	 * @param y1 The y coordinate of the first point
	 * @param x2 The x coordinate of the second point
	 * @param y2 The y coordinate of the second point
	 * @param thickness The thickness of the line
	 * @param stroke_color The color of the line
	 **/
	SVGLine(double x1, double y1, double x2, double y2, double thickness, const Color & stroke_color);

	// Methods

	/**
	 * @brief Implementation of the pure virtual function from SVGAnimatedElement
	 * @brief Outputs the SVGLine to the stream parameter 
	 * @param os The output stream
	 **/
	void print(std::ostream & os) const override;

private:
	// Members
	double m_x1;				// First point x coordinate
	double m_y1;				// First point y coordinate
	double m_x2;				// Second point x coordinate
	double m_y2;				// Second point y coordinate
	double m_thickness;			// Thickness of the line
	Color m_stroke_color;		// 
};