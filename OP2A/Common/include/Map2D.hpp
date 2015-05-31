/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * Map2D.hpp
 * 			-  
 *  
 */
#ifndef MAP2D_HPP_
#define MAP2D_HPP_


#include "Common/include/Exception_NoSuchValue.hpp"
#include "Common/include/Trio.hpp"
#include "Common/include/StringOps.hpp"


namespace OP2A{
namespace Common{


template <typename KEY1,	typename KEY2,	typename VALUE>
class Map2D
{
public:
	typedef typename std::vector< Trio<KEY1,KEY2,VALUE> >::iterator Iterator;

public:
	/*
	 * Constructor
	 */
	Map2D(size_t maxSize): m_sorted(false)
	{
		if(maxSize != 0)	m_vectorMap.reserve(maxSize);
	}


	/*
	 * Destructor
	 */
	~Map2D()
	{

	}



	/*
	 * =====================================================
	 * 	Internal Functions
	 * 	@author Minkwan KIM
	 * 	@version 1.0 /28/05/2015
	 * =====================================================
	 */

	// F01- Reserve
	// @brief	Reserve memory for the indicated size
	void reserve(const size_t* maxSize)
	{
		m_vectorMap.reserve(maxSize);
	}

	// F02- size
	// @brief	Return the size of this Map
	size_t size() const
	{
		return m_vectorMap.size();
	}

	// F03- begin()
	// @brief return the index of begin
	typename std::vector<Trio<KEY1, KEY2, VALUE> >::iterator begin()
	{
		return m_vectorMap.begin();
	}

	// F04- end()
	// @brief return the index of end
	typename std::vector<Trio<KEY1, KEY2, VALUE> >::iterator end()
	{
		return m_vectorMap.end();
	}

	// F05- insert
	// @beirf Insert value
	void insert(const KEY1& aKey1, const KEY2& aKey2, const VALUE& aValue)
	{
		  m_sorted = false;
		  m_vectorMap.push_back(Trio<KEY1, KEY2, VALUE>(aKey1, aKey2, aValue));
	}


	// F06- findBounds
	// @brief Find the bound of this map
	std::pair<typename std::vector<Trio<KEY1,KEY2,VALUE> >::iterator, typename std::vector<Trio<KEY1,KEY2,VALUE> >::iterator> findBounds(const KEY1& aKey1, const KEY2& aKey2)
	{
		if(m_vectorMap.empty())
		{
			valueNotFound(aKey1,aKey2);
		}

		if(!m_sorted)
		{
			sortKeys();
		}

		Iterator itr = std::lower_bound(m_vectorMap.begin(), m_vectorMap.end(), std::make_pair(aKey1,aKey2), LessThan());

		if((itr->first != aKey1) || (itr->second != aKey2))
		{
			valueNotFound(aKey1,aKey2);
		}

		Equal eq;
		std::pair<KEY1,KEY2> key(aKey1,aKey2);

		Iterator before = itr;
		for( ; itr != m_vectorMap.end(); ++itr)
		{
			if (!eq(*itr,key)) break;
		}

		return std::make_pair(before,itr);
	}


	// F07 Find
	// @brief find value
	VALUE find(const KEY1& aKey1, const KEY2& aKey2)
	{
		if(m_vectorMap.empty())
		{
			valueNotFound(aKey1,aKey2);
		}

		if(!m_sorted)
		{
			sortKeys();
		}

		Iterator itr = std::lower_bound(m_vectorMap.begin(), m_vectorMap.end(), std::make_pair(aKey1, aKey2), LessThan());

		if((itr->first != aKey1) || (itr->second != aKey2))
		{
			valueNotFound(aKey1,aKey2);
		}

		return itr->third;
	}

	VALUE find(const std::pair<KEY1,KEY2>& aKey)
	{
		return find(aKey.first, aKey.second);
	}


	// F08 sortKeys
	// @brief sorting keys
	void sortKeys()
	{
		std::sort(m_vectorMap.begin(),m_vectorMap.end(), LessThan());
		m_sorted = true;
	}


	// F09 is Present
	bool isPresent(const Trio<KEY1,KEY2,VALUE>& value)
	{
		op_assert(m_sorted);
		return std::binary_search(m_vectorMap.begin(), m_vectorMap.end(), value, LessThan());
	}

	// O01 - []
	// @return The value of this map.
	VALUE& operator[] (const OPuint i)
	{
		op_assert(i < size());
		return m_vectorMap[i].third;
	}


private:
	void valueNotFound(const KEY1& aKey1, const KEY2& aKey2)
	{
		std::string msg = "Map2D: KEYS not found: ";
		msg += StringOps::to_str(aKey1);
		msg += " ";
		msg += StringOps::to_str(aKey2);
		throw Common::ExceptionNoSuchValue (FromHere(),msg);
	}

private: // nested classes
	class LessThan
	{
	public:
		bool operator() (const Trio<KEY1,KEY2,VALUE>& p1, const Trio<KEY1,KEY2,VALUE>& p2) const
		{
			bool ret;
			if (p1.first == p2.first)
			{
				ret = (p1.second < p2.second);
			}
			else
			{
				ret = (p1.first < p2.first);
			}

			return ret;
		}

		bool operator() (const Trio<KEY1,KEY2,VALUE>& p1, const std::pair<KEY1,KEY2>& keys) const
		{
			bool ret;
			if (p1.first == keys.first)
			{
				ret = (p1.second < keys.second);
			}
			else
			{
				ret = (p1.first < keys.first);
			}
			return ret;
		}


		bool operator() (const std::pair<KEY1,KEY2>& keys, const Trio<KEY1,KEY2,VALUE>& p1) const
		{
			bool ret;
			if (p1.first == keys.first)
			{
				ret = (p1.second > keys.second);
			}
			else
			{
				ret = (p1.first > keys.first);
			}
			return ret;
		}
	}; // end class LessThan



	class Equal
	{
	public:
		bool operator() (const Trio<KEY1,KEY2,VALUE>& p1, const Trio<KEY1,KEY2,VALUE>& p2) const
		{
			return (p1.first == p2.first) && (p1.second == p2.second);
		}

		bool operator() (const Trio<KEY1,KEY2,VALUE>& p1, const std::pair<KEY1,KEY2>& keys) const
		{
			return (p1.first == keys.first) && (p1.second == keys.second);
		}

		bool operator() (const std::pair<KEY1,KEY2>& keys, const Trio<KEY1,KEY2,VALUE>& p1) const
		{
			return (p1.first == keys.first) && (p1.second == keys.second);
		}
	};

private:
	bool m_sorted;
	std::vector< Trio<KEY1, KEY2, VALUE> >  m_vectorMap;

};


}
}



#endif /* MAP2D_HPP_ */
