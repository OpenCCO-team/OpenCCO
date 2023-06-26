#pragma once

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include "svg_attributes.h"


namespace SVG
{

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               SVG::Animation                ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	/**
	 * @brief SVG animate element, inherits from Streamable
	 * @brief Convoluted class structure, but the logic is that an <animate> tag 
	 * @brief cannot have <animate> tags nested within
	 **/
	class Animation : public Streamable
	{
	public:
		// Constructors
		
		/**
		 * @brief Constructor, initializes all members with default constructor
		 **/
		Animation();

		/**
		 * @brief Constructor, initializes all members with parameters
		 **/
		Animation(const std::string & attribute_name,
			       const Timeline & timeline,
				   const int repeat_count = 0,
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
		 * @brief Outputs the Animation to the stream parameter 
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
////////////////////////             SVG::AnimatedElement            ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	/**
	 * @brief SVG animated element, inherits from Streamable
	 * @brief Abstract class, any SVG element that support animation will inherit from AnimatedElement
	 * @brief i.e elements that can have nested <animate> tags
	 **/
	class AnimatedElement : public Streamable
	{
	public:
		// Constructors

		/**
		 * @brief Constructor, initializes all members with default constructor
		 **/
		AnimatedElement();

		// Methods

		/**
		 * @brief Pure virtual method, output the AnimatedElement to the stream 
		 * @brief (called by operator<<(std::ostream & os, const Streamable & element) )
		 * @param os The output stream
		 **/
		virtual void print(std::ostream & os) const = 0;

		/**
		 * @brief Appends an SVG animation to m_animations
		 * @brief Checks if an animation for the same attribute exists or not, and by default will overwrite it
		 * @param animation The animation to append
		 **/
		void addAnimation(const Animation & animation);

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
		 * @returns a pointer to the Animation element if it exists, nullptr if it doesn't
		 **/
		Animation * getAnimation(const std::string & attribute_name);

		/**
		 * @brief Initialize an Animation for the AnimatedElement, for specified attribute and start & end values
		 * @param attribute_name The attribute to animate
		 * @param start_value The value of the attribute at the start of the timeline (keyTime 0.0)
		 * @param end_value The value of the attribute at the en of the timeline (keyTime 1.0)
		 * @param repeat_count How many time the animation will loop
		 * @param freeze Whether the animation freezes when it stops or if attribute resets
		 **/
		void initializeAnimation(const std::string & attribute_name,
								double start_value,
								double end_value,
								int duration,
								bool calc_mode_spline,
								int repeat_count = 0,
								bool freeze = false);

	private:
		std::vector<Animation> m_animations;
	};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  SVG::Svg                   ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	/**
	 * @brief SVG svg, inherits from AnimatedElement
	 * @brief Represents the <svg> ... </svg> tag
	 * @brief Also acts as a container of other svg elements
	 **/
	class Svg : public AnimatedElement
	{
	public:
		// Constructors

		/**
		 * @brief Constructor, initializes all members with default constructor
		 **/
		Svg();

		/**
		 * @brief Constructor, initializes all members with parameters
		 * @param x The x coordinate of the top left corner of the SVG
		 * @param y The y coordinate of the top left corner of the SVG
		 * @param width The width of the SVG
		 * @param height The height of the SVG
		 **/
		Svg(double x, double y, double width, double height);

		// Methods

		/**
		 * @brief Implementation of the pure virtual function from AnimatedElement
		 * @brief Outputs the Svg to the stream parameter 
		 * @param os The output stream
		 **/
		void print(std::ostream & os) const override;

		/**
		 * @brief Appends a shared pointer of an AnimatedElement to m_contained_elements
		 * @param ptr_element The shared pointer to append
		 **/
		void addElement(const std::shared_ptr<AnimatedElement> & ptr_element);

	private:
		// Members
		double m_x;					// Top left corner x coordinate
		double m_y;					// Top left corner y coordinate
		double m_width;				// SVG width
		double m_height;			// SVG height
		std::vector< std::shared_ptr<AnimatedElement> > m_contained_elements;	// svg elements within the svg tags
	};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  SVG::Rect                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	/**
	 * @brief SVG rectangle, inherits from AnimatedElement
	 * @brief Represents the <rect> ... </rect> tag
	 **/
	class Rect : public AnimatedElement
	{
	public:
		// Constructors

		/**
		 * @brief Constructor, initializes all members with default constructor
		 **/
		Rect();

		/**
		 * @brief Constructor, initializes all members with parameters
		 * @param x The x coordinate of the top left corner of the rectangle
		 * @param y The y coordinate of the top left corner of the rectangle
		 * @param width The width of the rectangle
		 * @param height The height of the rectangle
		 * @param fill_color The color of the rectangle
		 **/
		Rect(double x, double y, double width, double height, const Color & fill_color);

		// Methods

		/**
		 * @brief Implementation of the pure virtual function from AnimatedElement
		 * @brief Outputs the Rect to the stream parameter 
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
////////////////////////                  SVG::Line                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	/**
	 * @brief SVG line, inherits from AnimatedElement
	 * @brief Represents the <line> ... </line> tag
	 **/
	class Line : public AnimatedElement
	{
	public:
		// Constructors

		/**
		 * @brief Constructor, initializes all members with default constructor
		 **/
		Line();

		/**
		 * @brief Constructor, initializes all members with parameters
		 * @param x1 The x coordinate of the first point
		 * @param y1 The y coordinate of the first point
		 * @param x2 The x coordinate of the second point
		 * @param y2 The y coordinate of the second point
		 * @param thickness The thickness of the line
		 * @param stroke_color The color of the line
		 **/
		Line(double x1, double y1, double x2, double y2, double thickness, const Color & stroke_color);

		// Methods

		/**
		 * @brief Implementation of the pure virtual function from AnimatedElement
		 * @brief Outputs the Line to the stream parameter 
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
}