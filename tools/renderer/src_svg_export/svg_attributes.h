#pragma once

#include <vector>
#include <iostream>
#include <string>

#include "streamable.h"



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                    Color                    ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Color class to manage the RGB components
 * @brief Represents the color attribute in SVG
 **/
struct Color : public Streamable
{
	// Constructors

	/**
	 * @brief Constructor, initializes red, green and blue components to 0
	 **/
	Color();

	/**
	 * @brief Constructor, initializes red, green and blue components with parameters
	 * @param r The value of the red component
	 * @param g The value of the green component
	 * @param b The value of the blue component
	 **/
	Color(unsigned short r, unsigned short g, unsigned short b);

	// Methods

	/**
	 * @brief Implementation of the pure virtual function from Streamable
	 * @brief Outputs the Color to the stream parameter 
	 * @param os The output stream
	 **/
	void print(std::ostream & os) const override;

private:
	// Members
	unsigned short m_red;		// Red color component
	unsigned short m_green;		// Green color component
	unsigned short m_blue;		// Blue color component
};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   Timeline                  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Timeline class
 * @brief Intended for svg animate
 **/
struct Timeline : public Streamable
{
	// Constructors

	/**
	 * @brief Constructor, initializes members with default values
	 **/
	Timeline();

	/**
	 * @brief Constructor, initializes red, green and blue components with parameters
	 * @param r The value of the red component
	 * @param g The value of the green component
	 * @param b The value of the blue component
	 **/
	Timeline(int duration, bool calc_mode_spline);

	// Methods

	/**
	 * @brief Adds a key time and a value to the corresponding vector.
	 * @brief Keeps the vectors sorted by ascending keytime.
	 * @brief Replaces value if key_time is already in m_key_times
	 * @param key_time The key time (fraction of m_duration, between 0 and 1)
	 * @param value The value associated to the key time
	 **/
	void addKeyTime(const double key_time, const double value);

	/**
	 * @brief Adds a key time and a value to the corresponding vector.
	 * @brief Keeps the vectors sorted by ascending keytime.
	 * @brief Replaces value if key_time is already in m_key_times
	 * @param timestamp The timestamp of the key time (between 0 and m_duration)
	 * @param value The value associated to the key time
	 **/
	void addKeyTime(const int timestamp, const double value);

	/**
	 * @brief Implementation of the pure virtual function from Streamable
	 * @brief Outputs the Color to the stream parameter 
	 * @param os The output stream
	 **/
	void print(std::ostream & os) const override;

private:
	// Members
	int m_duration;						// Duration of one animation loop (in milliseconds)
	std::vector<double> m_key_times;	// Key animation times, each a fraction of the total duration (between 0 and 1)
	std::vector<double> m_values;		// values of the attribute at corresponding key times
	bool m_calc_mode_spline;			// Whether the animate attribute clacMode is equal to "spline" or not
};



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////               Other functions               ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Template function to output an SVG attribute to a stream
 * @brief Definition in header because it is a template function
 * @param os The output stream
 * @param attribute_name The name of the attribute
 * @param attribute_value The value associated to the attribute
 **/
template<class T>
void streamAttribute(std::ostream & os,
					const std::string& attribute_name,
					const T & attribute_value,
					const std::string & attribute_unit = "")
{
	os << attribute_name << "=\"" << attribute_value << attribute_unit << "\" ";
}