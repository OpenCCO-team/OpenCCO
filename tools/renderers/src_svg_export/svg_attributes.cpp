#include "svg_attributes.h"

#include <algorithm>



namespace SVG
{

/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                  SVG::Color                 ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	// Constructors


	Color::Color()
		: m_red(0), m_green(0), m_blue(0)
	{
		//
	}


	Color::Color(unsigned short r, unsigned short g, unsigned short b)
		: m_red(r), m_green(g), m_blue(b)
	{
		//
	}


	// Methods

	void Color::print(std::ostream & os) const
	{
		os << "rgb(" << m_red << "," << m_green << "," << m_blue << ")";
	}



/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 SVG::Timeline               ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


	// Constructors

	Timeline::Timeline()
	{
		//
	}


	Timeline::Timeline(int duration, bool calc_mode_spline)
		: m_duration(duration), m_calc_mode_spline(calc_mode_spline)
	{
		//
	}


	// Methods

	void Timeline::print(std::ostream & os) const
	{
		// duration="<m_duration>ms"
		streamAttribute(os, "dur", m_duration, "ms");

		// values="<m_values[0]; ..."
		os << "values=\"";

		std::vector<double>::const_iterator itval = m_values.begin();
		if(itval != m_values.end())
		{
			os << *(itval++);
		}

		while(itval != m_values.end())
		{
			os << "; " << *(itval++);
		}

		os << "\" ";

		// keyTimes="<m_key_times[0]; ..."
		os << "keyTimes=\"";

		std::vector<double>::const_iterator itkeytime = m_key_times.begin();
		if(itkeytime != m_key_times.end())
		{
			os << *(itkeytime++);
		}

		while(itkeytime != m_key_times.end())
		{
			os << "; " << *(itkeytime++);
		}

		os << "\" ";

		// optional spline calc
		if(m_calc_mode_spline && m_values.size() >= 2)	// at least 2 values to use splines
		{
			// calcMode="spline"
			streamAttribute(os, "calcMode", "spline");

			std::string spline = "0.5 0 0.5 1";			// soft curve from 0 to 1

			// there's one less spline than key time or value (they're used to join them)
			// keySplines="<spline>; ..."

			os << "keySplines=\"" << spline;	// first iteration, not prefixed by "; "
			for(int i = 1; i < m_values.size() - 1; i++)	// starting at i=1 because we already did first
			{
				os << "; " << spline;
			}
			os << "\" ";
		}
	}


	void Timeline::addKeyTime(const double key_time, const double value)
	{
		// only add the key if within [0; 1]
		if(key_time >= 0.0 && key_time <= 1.0)
		{
			// lamba to test if >= value
			auto supeq_keyTime = [&key_time] (double x) { return x >= key_time; };

			std::vector<double>::iterator it = std::find_if(m_key_times.begin(), m_key_times.end(), supeq_keyTime);

			if(it == m_key_times.end())		// key_time highest key time
			{
				m_key_times.push_back(key_time);
				m_values.push_back(value);
			}
			else							// key_time must be inserted (or it's value edited if already in)
			{
				int index = it - m_key_times.begin();
				// it points to the first element >= key_time
				if(*it != key_time)			// the key time at the it position is > key_time
				{
					m_key_times.insert(m_key_times.begin() + index, key_time);
					m_values.insert(m_values.begin() + index, value);
				}
				else						// the key_time already exists; we update the associated value
				{
					std::vector<double>::iterator itval = m_values.begin() + index;

					*itval = value;
				}
			}
		}
		else
		{
			// do nothing, key_time out of bounds
			// TODO : error handling
		}
	}

	void Timeline::addKeyTime(const int timestamp, const double value)
	{
		// valid time is anywhere between 0 and m_duration
		if(m_duration != 0 && timestamp >= 0 && timestamp <= m_duration)
		{
			double key_time = ((double) timestamp) / (double) m_duration;
			addKeyTime(key_time, value);
		}
	}
}
