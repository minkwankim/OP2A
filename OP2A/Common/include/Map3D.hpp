/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2015
 *      			Author: Minkwan Kim
 *
 * Map3D.hpp
 * 			-  
 *  
 */
#ifndef MAP3D_HPP_
#define MAP3D_HPP_


#include "Common/include/OP2A.hpp"
#include "Common/include/Exception_NoSuchValue.hpp"
#include "Common/include/StringOps.hpp"
#include "Common/include/Quartet.hpp"



namespace OP2A{
namespace Common{

template <typename KEY1, typename KEY2, typename KEY3, typename VALUE>
class Map3D
{
public:
	typedef typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator MapIterator;

public:
	Map3D(size_t maxSize):m_sorted(false)
	{
		if(maxSize != 0)
		{
			m_vectorMap.reserve(maxSize);
		}
	}

	/// Default destructor
	~Map3D()
	{

	}




	/*
	 * =======================================
	 * Internal Functions
	 * =======================================
	 */
	// F01 - reserve
	void reserve(const size_t& maxSize)
	{
		m_vectorMap.reserve(maxSize);
	}


	// F02 - clear
	// Clear the content of the map
	void clear()
	{
		std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >().swap(m_vectorMap);
	}


	// F03 - insert
	/// Insert a value with corresponding keys
	void insert(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3, const VALUE& aValue)
	{
		m_sorted = false;
		m_vectorMap.push_back(Quartet<KEY1, KEY2, KEY3, VALUE>(aKey1, aKey2, aKey3, aValue));
	}

	void insert(const Trio<KEY1,KEY2,KEY3>& aKey, const VALUE& aValue)
	{
		insert(aKey.first, aKey.second, aKey.third, aValue);
	}


	// F04 - findBounds
	/// Find the upper and lower pairs of KEY'S and VALUE'S bounding the supplied KEY.
	std::pair<typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator, typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator>
	findBounds(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
	{
		if(m_vectorMap.empty())
		{
			valueNotFound(aKey1,aKey2,aKey3);
		}

		if(!m_sorted)
		{
			sortKeys();
		}

		MapIterator itr = std::lower_bound(m_vectorMap.begin(), m_vectorMap.end(), make_Trio(aKey1,aKey2,aKey3), LessThan());

		if((itr->first != aKey1) || (itr->second != aKey2) || (itr->third != aKey3))
		{
			valueNotFound(aKey1,aKey2,aKey3);
		}

		Equal eq;
		Trio<KEY1,KEY2,KEY3> key(aKey1,aKey2,aKey3);
		MapIterator before = itr;

		for( ; itr != m_vectorMap.end(); ++itr)
		{
			if (!eq(*itr,key)) break;
		}

		return std::make_pair(before,itr);
	}


	// F05 - find
	// Find VALUE with the given KEY
	VALUE find(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
	{
		if(m_vectorMap.empty())
		{
			valueNotFound(aKey1,aKey2,aKey3);
		}

		if(!m_sorted)
		{
			sortKeys();
		}

		MapIterator itr = std::lower_bound(m_vectorMap.begin(), m_vectorMap.end(), make_Trio(aKey1,aKey2,aKey3), LessThan());

		if((itr->first != aKey1) || (itr->second != aKey2) || (itr->third != aKey3))
		{
			valueNotFound(aKey1,aKey2,aKey3);
		}

		return itr->fourth;
	}


	std::pair<MapIterator, MapIterator> find(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3, bool& isFound)
	{
		if(!m_sorted)
		{
			sortKeys();
		}

		std::pair<MapIterator,MapIterator> result;

		if (m_vectorMap.empty())
		{
			isFound = false;
			return result;
		}

		result=std::equal_range(m_vectorMap.begin(), m_vectorMap.end(), make_Trio(aKey1,aKey2,aKey3), LessThan());

		isFound = !((result.first->first != aKey1) || (result.first->second != aKey2) || (result.first->third != aKey3));

		return result;
	}

	VALUE find(const Trio<KEY1,KEY2,KEY3>& aKey)
	{
		return find(aKey.first, aKey.second, aKey.third);
	}


	// F06 - size
	// Get the number of pairs already inserted
	size_t size() const
	{
		return m_vectorMap.size();
	}


	// F07 - print
	// Prints all the elements of the map to the ostream
	std::ostream& print(std::ostream& out)
	{
		for(OPuint i = 0; i < m_vectorMap.size(); ++i)
		{
			out << m_vectorMap[i] << "\n";
		}
		return out;
	}


	// F08 - []
	// Overloading of the operator"[]" for assignment
	VALUE& operator[] (const OPuint i)
	{
		op_assert(i < size());
		return m_vectorMap[i].fourth;
	}


	// F09 - getEntry
	// Gets the i entry
	Quartet<KEY1,KEY2,KEY3,VALUE>& getEntry(const OPuint i)
	{
		op_assert(i < size());
		return m_vectorMap[i];
	}


	// F10 - sortKeys
	// Sort all the pairs in the map by key
	void sortKeys()
	{
		std::sort(m_vectorMap.begin(),m_vectorMap.end(), LessThan());
		m_sorted = true;
	}


	// F11 - exists
	// Checks if an element with a certain key has been inserted
	bool exists (const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
	{
		if(!m_sorted)
		{
			sortKeys();
		}

		MapIterator itr = std::find_if ( m_vectorMap.begin(), m_vectorMap.end(), EqualComp ( make_Trio(aKey1,aKey2,aKey3 ) ) );
		return itr != m_vectorMap.end();
	}


	// F12 - begin
	typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator begin()
	{
		return m_vectorMap.begin();
	}


	// F13 - end
	typename std::vector< Quartet<KEY1,KEY2,KEY3,VALUE> >::iterator end()
	{
		return m_vectorMap.end();
	}


private:
	/// Error reporting in case of the key not being present
	void valueNotFound(const KEY1& aKey1, const KEY2& aKey2, const KEY3& aKey3)
	{
		std::string msg = "CFMap3D: KEYS not found: ";
		msg += Common::StringOps::to_str(aKey1);
		msg += " ";
		msg += Common::StringOps::to_str(aKey2);
		msg += " ";
		msg += Common::StringOps::to_str(aKey3);
		throw Common::ExceptionNoSuchValue (FromHere(),msg);
	}


	/*
	 * ==========================
	 * Nested Class
	 * ==========================
	 */
private:
	class LessThan
	{
	public:
		bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1, const Quartet<KEY1,KEY2,KEY3,VALUE>& p2) const
		{
			bool ret;
			if (p1.first == p2.first)
			{
				if (p1.second == p2.second)
				{
					ret = (p1.third < p2.third);
				}
				else
				{
					ret = (p1.second < p2.second);
				}
			}
			else
			{
				ret = (p1.first < p2.first);
			}

			return ret;
		}

		bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1, const Trio<KEY1,KEY2,KEY3>& keys) const
		{
			bool ret;
			if (p1.first == keys.first)
			{
				if (p1.second == keys.second)
				{
					ret = (p1.third < keys.third);
				}
				else
				{
					ret = (p1.second < keys.second);
				}
			}
			else
			{
				ret = (p1.first < keys.first);
			}

			return ret;
		}


		bool operator() (const Trio<KEY1,KEY2,KEY3>& keys, const Quartet<KEY1,KEY2,KEY3,VALUE>& p1) const
		{
			bool ret;
			if (p1.first == keys.first)
			{
				if (p1.second == keys.second)
				{
					ret = (p1.third > keys.third);
				}
				else
				{
					ret = (p1.second > keys.second);
				}
			}
			else
			{
				ret = (p1.first > keys.first);
			}

			return ret;
		}
	};


	class Equal
	{
	public:
		bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1, const Quartet<KEY1,KEY2,KEY3,VALUE>& p2) const
		{
			return (p1.first == p2.first) && (p1.second == p2.second) && (p1.third == p2.third);
		}

		bool operator() (const Quartet<KEY1,KEY2,KEY3,VALUE>& p1, const Trio<KEY1,KEY2,KEY3>& keys) const
		{
			return operator()(keys,p1);
		}

		bool operator() (const Trio<KEY1,KEY2,KEY3>& keys, const Quartet<KEY1,KEY2,KEY3,VALUE>& p1) const
		{
			return (p1.first == keys.first) && (p1.second == keys.second) && (p1.third == keys.third);
		}

	}; // end class Equal

	struct EqualComp
	{
		const Trio<KEY1,KEY2,KEY3>& m_key;
		EqualComp (const Trio<KEY1,KEY2,KEY3>& aKey) : m_key (aKey) {}

		bool operator() ( const Quartet<KEY1,KEY2,KEY3,VALUE>& obj )
		{
			return (m_key.first == obj.first) && (m_key.second == obj.second) && (m_key.third == obj.third);
		};
	};


//////////////////////////////////////////////////////////////////////////////





private: //data
  bool m_sorted;
  std::vector< Quartet<KEY1, KEY2, KEY3, VALUE> > m_vectorMap;


};



}
}


#endif /* MAP3D_HPP_ */
