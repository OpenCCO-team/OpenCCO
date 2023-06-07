#pragma once

#include <iostream>



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  Streamable                 ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Abstract class, any derived class can be streamed with operator<<
 **/
class Streamable
{
public:
	// Constructors
	
	/**
	 * @brief Constructor
	 **/
	Streamable();

	// Methods

	/**
	 * @brief Pure virtual method, output the Streamable to the stream 
	 * @brief (called by operator<<(std::ostream & os, const Streamable & element) )
	 * @param os The output stream
	 **/
	virtual void print(std::ostream & os) const = 0;

	// Operator overloads

	/**
	 * @brief operator<< overload for displaying and file writing purpose
	 * @param os The output stream
	 * @param color The SVG element to output
	 * @returns The output stream
	 **/
	friend std::ostream& operator<<(std::ostream & os, const Streamable & s);
};