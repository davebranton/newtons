/**
	argstream.h
	Purpose: Parse the command line. It's a mod of original argstream developed
	by Xavier Décoret. This mod adds support for muti-byte argument string and
	wide-char argument string.
	Copyright (C) 2004 Xavier Décoret <Xavier.Decoret@imag.fr>
	Portions Copyright (C) 2015 Levski Weng <levskiweng@gmail.com>

	@author Xavier Décoret
	@author Levski Weng
	@version 1.0

	@copyright argsteam is a free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	@copyright argstream is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	@copyright You should have received a copy of the GNU General Public License
	along with argstream; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef ARGSTREAM_H
#define ARGSTREAM_H

#include <string>
#include <list>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <locale>
#include <codecvt>
#include <type_traits>
#include <cassert>
#include <memory>

namespace argstream
{
	//--------------------------------------------------------------------------

	typedef enum
    {
		PARSED_OK = 0,
		PARSED_ERR_HELP_REQUESTED,
		PARSED_ERR_UNUSED_PARAMETER,
		PARSED_ERR_OTHER
	} RESULT_OF_PARSE;

	/**
       Main class to store the argument string.
	*/
	template<typename CHARTYPE>
    class argstream;

	/**
       Get the help description.
	*/
	template<typename CHARTYPE, typename T>
    struct description_policy;

	/**
       The option holder which store the specified option.
	*/
	template <typename CHARTYPE>
    class OptionHolder;

	/**
       The value holder which stores the value of the speicified argument.
	*/
	template<typename CHARTYPE, typename T>
    class ValueHolder;

	/*
	template<typename CHARTYPE, typename T, typename O>
    class ValuesHolder;
	*/

	/**
       Convert UTF-8 string to UTF-16 and vice versa.
	*/
	template<typename CHARTYPE>
    struct TSTR;

	/**
       Template string stream
	*/
	template<typename CHARTYPE>
    struct TSTRSTREAM;

	/**
       Store the copyright information.
	*/
	template <typename CHARTYPE>
    class CopyrightHolder;

	/**
       Store the command usage examples.
	*/
	template <typename CHARTYPE>
    class ExampleHolder;

	/**
       Parse the command line and store the specified parameter value.
	*/
	template <typename CHARTYPE, typename T>
	inline ValueHolder<CHARTYPE, T>
	parameter(
              CHARTYPE s,
              const CHARTYPE* l,
              T& b,
              const CHARTYPE* desc,
              bool mandatory = true
             );

	/* Disable it temporarily
	template<typename CHARTYPE, typename T, typename O>
	inline ValuesHolder<CHARTYPE, T, O>
	values(const O& o, const CHARTYPE* desc, int len=-1);
	*/

	/**
		Generate the option.

		@param s Short parameter name.
		@param l Long parameter name.
		@param b Whether the parameter appears or not.
		@param desc The description of the parameter.

		@return The option holder.
	*/
	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>
	option(
           CHARTYPE s,
           const CHARTYPE* l,
           bool& b,
           const CHARTYPE* desc
          );

	/**
		Generate the help option.
	*/
	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>
	help();

	/**
		Generate the examples.

		@param cmdline Example command line.
		@param desc The description of the command line.

		@return The object that holds the example.
	*/
	template<typename CHARTYPE>
	inline ExampleHolder<CHARTYPE>
	example(
            const CHARTYPE* cmdline,
            const CHARTYPE* desc
           );

	/**
		Generate the copyright information which will appear on the top of the help message.

		@param copyright Copyright message.

		@return The object that holds the copyright information.
	*/
	template<typename CHARTYPE>
	inline CopyrightHolder<CHARTYPE>
	copyright(const CHARTYPE* copyright);

	/**
		Parse the "option - value" parameter.

		@param s Reference to the argstream object which is going to be parsed.
		@param v Reference to the value holder which stores the value.

		@return Reference to the parsed argstream object.
	*/
	template<typename CHARTYPE, typename T>
	inline argstream<CHARTYPE>&
	operator >>(argstream<CHARTYPE>& s, ValueHolder<CHARTYPE, T> const& v);

	/**
		Parse the "option - value1 value2 value3" parameters.

		@param s Reference to the argstream object which is going to be parsed.
		@param v Reference to the values holder which stores the values.

		@return Reference to the parsed argstream object.
	*/
	/*
	template<typename CHARTYPE, typename T, typename O>
	inline argstream<CHARTYPE>&
	operator >>(argstream<CHARTYPE>& s, ValuesHolder<CHARTYPE, T, O> const& v);
	*/

	/**
		Parse the "option" parameters.

		@param s Reference to the argstream object which is going to be parsed.
		@param v Reference to the option holder which stores the existence of the option.

		@return Reference to the parsed argstream object.
	*/
	template<typename CHARTYPE>
	inline argstream<CHARTYPE>&
	operator >>(argstream<CHARTYPE>& s, OptionHolder<CHARTYPE> const& v);

	/**
		Add examples section to the help output.

		@param s Reference to the argstream object which is going to add the example.
		@param e Reference to the object which holds the example.

		@return Reference to the argstream object which has added the example.
	*/
	template<typename CHARTYPE>
	inline argstream<CHARTYPE>&
    operator >>(argstream<CHARTYPE>& s, ExampleHolder<CHARTYPE> const& e);

	/**
		Add copyright section to the help output.

		@param s Reference to the argstream object which is going to add the copyright information.
		@param c Reference to the object which holds the copyright information.

		@return Reference to the argstream object which has added the copyright information.
	*/
	template<typename CHARTYPE>
	inline argstream<CHARTYPE>&
    operator >>(argstream<CHARTYPE>& s, CopyrightHolder<CHARTYPE> const& c);
	//--------------------------------------------------------------------------

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of TSTR<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE>
	struct TSTR
	{
		typedef std::basic_string<CHARTYPE, std::char_traits<CHARTYPE>, std::allocator<CHARTYPE>> type;

		static inline typename TSTR::type ToString(const char* utf8_str);
		static inline typename TSTR::type ToString(const wchar_t* wstr);
		static inline typename TSTR::type ToString(char c);
		static inline typename TSTR::type ToString(wchar_t wc);
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of TSTR<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<>
	std::wstring TSTR<wchar_t>::ToString(const char* utf8_str)
	{
		std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
		return conv.from_bytes(utf8_str);
	}

	template<>
	std::wstring TSTR<wchar_t>::ToString(const wchar_t* wstr)
	{
		return std::wstring(wstr);
	}

	template<>
	std::wstring TSTR<wchar_t>::ToString(char c)
	{
		std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
		return conv.from_bytes(c);
	}

	template<>
	std::wstring TSTR<wchar_t>::ToString(wchar_t wc)
	{
		wchar_t buf[] = {wc, 0};
		return std::wstring(buf);
	}

	template<>
	std::string TSTR<char>::ToString(const char* utf8_str)
	{
		return std::string(utf8_str);
	}

	template<>
	std::string TSTR<char>::ToString(const wchar_t* wstr)
	{
		std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
		return conv.to_bytes(wstr);
	}

	template<>
	std::string TSTR<char>::ToString(char c)
	{
		char buf[] = {c, 0};
		return std::string(buf);
	}

	template<>
	std::string TSTR<char>::ToString(wchar_t wc)
	{
		std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
		return conv.to_bytes(wc);
	}

	/**
		@todo Expand this struct to hold file string and etc.
	*/
	template<typename CHARTYPE>
	struct TSTRSTREAM
	{
		typedef std::basic_istringstream<CHARTYPE, std::char_traits<CHARTYPE>,std::allocator<CHARTYPE>> I;
		typedef std::basic_ostringstream<CHARTYPE, std::char_traits<CHARTYPE>, std::allocator<CHARTYPE>> O;
		typedef std::basic_ostream<CHARTYPE, std::char_traits<CHARTYPE> > COUT;
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of ValueHolder<CHARTYPE, T>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE, typename T>
	class ValueHolder
	{
	public:
		ValueHolder(CHARTYPE s,
			const CHARTYPE* l,
			T& b,
			const CHARTYPE* desc,
			bool mandatory);
		ValueHolder(const CHARTYPE* l,
			T& b,
			const CHARTYPE* desc,
			bool mandatory);
		ValueHolder(CHARTYPE s,
			T& b,
			const CHARTYPE* desc,
			bool mandatory);

		template<typename C, typename T2>
		friend argstream<C>& operator>>(argstream<C>& s, ValueHolder<C, T2> const& v);

		friend struct description_policy<CHARTYPE, T>;

		typename TSTR<CHARTYPE>::type name() const;
		typename TSTR<CHARTYPE>::type description() const;
	private:
		typename TSTR<CHARTYPE>::type shortName_;
		typename TSTR<CHARTYPE>::type longName_;
		T* value_;
		T initialValue_;
		typename TSTR<CHARTYPE>::type description_;
		bool mandatory_;
	};
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of ValueHodler<CHARTYPE, T>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE, typename T>
	ValueHolder<CHARTYPE, T>::ValueHolder(
		CHARTYPE s,
		const CHARTYPE* l,
		T& v,
		const CHARTYPE* desc,
		bool mandatory)
		:  shortName_(1,s),
		longName_(l),
		value_(&v),
		initialValue_(v),
		description_(desc),
		mandatory_(mandatory)
	{
	}
	template<typename CHARTYPE, typename T>
	ValueHolder<CHARTYPE, T>::ValueHolder(
		const CHARTYPE* l,
		T& v,
		const CHARTYPE* desc,
		bool mandatory)
		:  longName_(l),
		value_(&v),
		initialValue_(v),
		description_(desc),
		mandatory_(mandatory)
	{
	}
	template<typename CHARTYPE, typename T>
	ValueHolder<CHARTYPE, T>::ValueHolder(CHARTYPE s,
		T& v,
		const CHARTYPE* desc,
		bool mandatory)
		:  shortName_(1,s),
		value_(&v),
		initialValue_(v),
		description_(desc),
		mandatory_(mandatory)
	{
	}
	template<typename CHARTYPE, typename T>
	typename TSTR<CHARTYPE>::type ValueHolder<CHARTYPE, T>::name() const
	{
		typename TSTRSTREAM<CHARTYPE>::O os;
		if (!shortName_.empty()) os<<'-'<<shortName_;
		if (!longName_.empty()) {
			if (!shortName_.empty()) os<<'/';
			os<<TSTR<CHARTYPE>::ToString("--")<<longName_;
		}
		return os.str();
	}

	template<typename CHARTYPE, typename T>
	inline typename TSTR<CHARTYPE>::type ValueHolder<CHARTYPE, T>::description() const
	{
		return description_;
	}


	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of CopyrightHolder<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template <typename CHARTYPE>
	class CopyrightHolder
	{
	public:
		inline CopyrightHolder(const CHARTYPE* copyright);
		inline typename TSTR<CHARTYPE>::type copyright() const;

		template<typename C>
		friend argstream<C>& operator>>(argstream<C>& s, CopyrightHolder<C> const& e);
	private:
		typename TSTR<CHARTYPE>::type copyright_;
	};
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of CopyrightHolder<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE>
	inline CopyrightHolder<CHARTYPE>::CopyrightHolder(const CHARTYPE* copyright)
		:copyright_(copyright)
	{
	}

	template<typename CHARTYPE>
	inline typename TSTR<CHARTYPE>::type CopyrightHolder<CHARTYPE>::copyright() const
	{
		return copyright_;
	}

	template<typename CHARTYPE>
	inline argstream<CHARTYPE>&
    operator >>(argstream<CHARTYPE>& s, CopyrightHolder<CHARTYPE> const& v)
	{
		s.copyright_ = v.copyright();
		return s;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of ExampleHolder<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template <typename CHARTYPE>
	class ExampleHolder
	{
	public:
		ExampleHolder(const CHARTYPE* cmdline, const CHARTYPE* desc);
		typename TSTR<CHARTYPE>::type cmdline() const;
		typename TSTR<CHARTYPE>::type description() const;
		template<typename C>
		friend argstream<C>& operator>>(argstream<C>& s, ExampleHolder<C> const& e);
	private:
		typename TSTR<CHARTYPE>::type cmdline_;
		typename TSTR<CHARTYPE>::type description_;
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of ExampleHolder<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE>
	ExampleHolder<CHARTYPE>::ExampleHolder(const CHARTYPE* cmdline, const CHARTYPE* desc)
		:cmdline_(cmdline), description_(desc)
	{
	}

	template<typename CHARTYPE>
	typename TSTR<CHARTYPE>::type
    ExampleHolder<CHARTYPE>::cmdline() const
	{
		return cmdline_;
	}

	template<typename CHARTYPE>
	typename TSTR<CHARTYPE>::type ExampleHolder<CHARTYPE>::description() const
	{
		return description_;
	}

	template<typename CHARTYPE>
	inline argstream<CHARTYPE>&
    operator >>(argstream<CHARTYPE>& s, ExampleHolder<CHARTYPE> const& v)
	{
		s.argExamples_.push_back(argstream<CHARTYPE>::example_entry(v.cmdline(), v.description()));
		return s;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of OptionHolder<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template <typename CHARTYPE>
	class OptionHolder
	{
	public:
		inline OptionHolder(
			CHARTYPE s,
			const CHARTYPE* l,
			bool& b,
			const CHARTYPE* desc);
		inline OptionHolder(
			const CHARTYPE* l,
			bool& b,
			const CHARTYPE* desc);
		inline OptionHolder(
			CHARTYPE s,
			bool& b,
			const CHARTYPE* desc);
		inline OptionHolder(CHARTYPE s,
			const CHARTYPE* l,
			const CHARTYPE* desc);

		inline typename TSTR<CHARTYPE>::type name() const;
		inline typename TSTR<CHARTYPE>::type description() const;

		template<typename C>
		friend argstream<C>& operator>>(argstream<C>& s, OptionHolder<C> const& v);

		friend OptionHolder<CHARTYPE> help<CHARTYPE>();
	private:
		typename TSTR<CHARTYPE>::type shortName_;
		typename TSTR<CHARTYPE>::type longName_;
		bool* value_;
		typename TSTR<CHARTYPE>::type description_;
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of OptionHolder<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>::OptionHolder(
		CHARTYPE s,
		const CHARTYPE* l,
		bool& b,
		const CHARTYPE* desc)
		: shortName_(1,s),
		longName_(l),
		value_(&b),
		description_(desc)
	{
	}

	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>::OptionHolder(
		const CHARTYPE* l,
		bool& b,
		const CHARTYPE* desc)
		: longName_(l),
		value_(&b),
		description_(desc)
	{
	}

	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>::OptionHolder(
		CHARTYPE s,
		bool& b,
		const CHARTYPE* desc)
		: shortName_(1,s),
		value_(&b),
		description_(desc)
	{
	}

	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>::OptionHolder(
		CHARTYPE s,
		const CHARTYPE* l,
		const CHARTYPE* desc)
		: shortName_(1,s),
		longName_(l),
		value_(NULL),
		description_(desc)
	{
	}

	template<typename CHARTYPE>
	inline typename TSTR<CHARTYPE>::type OptionHolder<CHARTYPE>::name() const
	{
		typename TSTRSTREAM<CHARTYPE>::O os;
		if (!shortName_.empty()) os<<'-'<<shortName_;
		if (!longName_.empty())
		{
			if (!shortName_.empty()) os<<'/';
			os<<TSTR<CHARTYPE>::ToString("--")<<longName_;
		}
		return os.str();
	}

	template<typename CHARTYPE>
	inline typename TSTR<CHARTYPE>::type OptionHolder<CHARTYPE>::description() const
	{
		return description_;
	}

	/* ValusesHolder is not used currently, disable it temporarily
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of ValuesHolder<CHARTYPE, T, O>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE, typename T, typename O>
	class ValuesHolder
	{
	public:
		ValuesHolder(const O& o,
			const CHARTYPE* desc,
			int len);

		template<typename C, typename T2, typename O2>
		friend argstream<C>& operator>>(argstream<C>& s, ValuesHolder<C, T2, O2> const& v);
		typename TSTR<CHARTYPE>::type name() const;
		typename TSTR<CHARTYPE>::type description() const;
		typedef T value_type;
	private:
		mutable O   value_;
		typename TSTR<CHARTYPE>::type description_;
		int len_;
		CHARTYPE letter_;
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of ValuesHolder<CHARTYPE, T, O>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE, typename T, typename O>
	ValuesHolder<CHARTYPE, T, O>::ValuesHolder(
		const O& o,
		const CHARTYPE* desc,
		int len)
		: value_(o),
		description_(desc),
		len_(len)
	{
		letter_ = argstream::uniqueLetter();
	}

	template <typename CHARTYPE, typename T, typename O>
	typename TSTR<CHARTYPE>::type
	ValuesHolder<CHARTYPE, T, O>::name() const
	{
		TSTRSTREAM<CHARTYPE>::O os;
		os<< letter_ <<TSTR<CHARTYPE>::ToString("i");
		return os.str();
	}

	template <typename CHARTYPE, typename T, typename O>
	typename TSTR<CHARTYPE>::type
	ValuesHolder<CHARTYPE, T, O>::description() const
	{
		return description_;
	}

	template<typename CHARTYPE, typename T, typename O>
	argstream<CHARTYPE>& operator >>(argstream<CHARTYPE>& s,ValuesHolder<CHARTYPE,T,O> const& v)
	{
		s.argHelps_.push_back(argstream::help_entry(v.name(),v.description()));
		{
			STRSTREAM<CHARTYPE>::O os;
			os<<' '<<v.letter_<<'1';
			switch (v.len_)
			{
			case -1:
				os<< TSTR<CHARTYPE>::ToString("...");
				break;
			case 1:
				break;
			default:
				os<< TSTR<CHARTYPE>::ToString("...") << v.letter_ << v.len_;
				break;
			}
			s.cmdLine_ += os.str();
		}
		std::list<TSTR<CHARTYPE>::type>::iterator first = s.values_.begin();
		// We add to the iterator as much values as we can, limited to the length
		// specified (if different of -1)
		int n = v.len_ != -1?v.len_:s.values_.size();
		while (first != s.values_.end() && n-->0)
		{
			// Read the value from the string *first
			ValueParser<CHARTYPE, T> p;
			*(v.value_++) = p(*first );
			s.argHelps_.push_back(argstream::help_entry(v.name(),v.description()));
			// The value we just removed was maybe "remembered" by an option so we
			// remove it now.
			for (std::map<typename TSTR<CHARTYPE>::type, typename argstream<CHARTYPE>::value_iterator>::iterator
				jter = s.options_.begin();jter != s.options_.end();++jter)
			{
				if (jter->second == first)
				{
					jter->second = s.values_.end();
				}
			}
			++first;
		}
		// Check if we have enough values
		if (n != 0)
		{
			s.isOk_ = false;
			typename TSTRSTREAM<CHARTYPE>::O os;
			os<<TSTR<CHARTYPE>::ToString("Expecting ")<<v.len_<<TSTR<CHARTYPE>::ToString(" values");
			s.errors_.push_back(os.str());
		}
		// Erase the values parsed
		s.values_.erase(s.values_.begin(),first);
		return s;
	}
	*/
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface and implementation of ValueParser<CHARTYPE, T>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE, typename T>
	class ValueParser
	{
	public:
		inline T operator ()(const typename TSTR<CHARTYPE>::type& s) const
		{
			typename TSTRSTREAM<CHARTYPE>::I is(s);
			T t;
			is>>t;
			return t;
		}
	};

	template<typename CHARTYPE>
	class ValueParser<CHARTYPE, typename TSTR<CHARTYPE>::type>
	{
	public:
		inline typename TSTR<CHARTYPE>::type
		operator ()(const typename TSTR<CHARTYPE>::type& s) const
		{
			return s;
		}
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Interface of argstream<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE>
	class argstream
	{
	public:
		inline argstream<CHARTYPE>(int argc,CHARTYPE const* const argv[]);
		inline argstream<CHARTYPE>(const CHARTYPE* c);

		template<typename C, typename T>
		friend argstream<C>& operator>>(
			argstream<C>& s,
			ValueHolder<C, T> const& v);
	
		template<typename C>
		friend argstream<C>& operator>>(
			argstream<C>& s,
			OptionHolder<C> const& v);
		/*
		template<typename T, typename O>
		friend argstream<CHARTYPE>& operator>>(
			argstream<CHARTYPE>& s,
			ValuesHolder<CHARTYPE, T,O> const& v);
		*/
		
		template<typename C>
		friend argstream<C>& operator>>(
			argstream<C>& s,
			ExampleHolder<C> const& v);

		template<typename C>
		friend argstream<C>& operator>>(
			argstream<C>& s,
			CopyrightHolder<C> const& v);

		inline bool helpRequested() const;
		inline bool isOk() const;
		inline typename TSTR<CHARTYPE>::type errorLog() const;
		inline typename TSTR<CHARTYPE>::type usage() const;
		inline RESULT_OF_PARSE defaultErrorHandling(bool ignoreUnused=false) const;
	protected:
		void parse(int argc, CHARTYPE const* const argv[]);
	private:
		typedef CHARTYPE* PCHARTYPE;
		typedef typename std::list<typename TSTR<CHARTYPE>::type>::iterator value_iterator;
		typedef typename std::pair<typename TSTR<CHARTYPE>::type, typename TSTR<CHARTYPE>::type> help_entry;
		typedef typename std::pair<typename TSTR<CHARTYPE>::type, typename TSTR<CHARTYPE>::type> example_entry;
		typename TSTR<CHARTYPE>::type progName_;
		typename TSTR<CHARTYPE>::type cmdLine_;
		typename TSTR<CHARTYPE>::type copyright_;
		std::map<typename TSTR<CHARTYPE>::type, value_iterator> options_;
		std::list<typename TSTR<CHARTYPE>::type> values_;
		bool minusActive_;
		bool isOk_;
		std::unique_ptr<PCHARTYPE> argv_from_cmdline_;
		std::deque<std::pair<typename TSTR<CHARTYPE>::type, typename TSTR<CHARTYPE>::type>> argHelps_;
		std::deque<std::pair<typename TSTR<CHARTYPE>::type, typename TSTR<CHARTYPE>::type>> argExamples_;
		std::deque<typename TSTR<CHARTYPE>::type> errors_;
		bool helpRequested_;
	};

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of argstream<CHARTYPE>
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template<typename CHARTYPE>
	inline argstream<CHARTYPE>::argstream(int argc, CHARTYPE const* const argv[])
		: progName_(),
		minusActive_(true),
		isOk_(true),
		argv_from_cmdline_(nullptr),
		helpRequested_(false)
	{
		typename TSTR<CHARTYPE>::type argv0(argv[0]);
		size_t found = argv0.find_last_of(TSTR<CHARTYPE>::ToString("/\\"));
		progName_ = argv0.substr(found+1);
		parse(argc,argv);
	}

	template<typename CHARTYPE>
	inline argstream<CHARTYPE>::argstream(const CHARTYPE* c)
		: progName_(TSTR<CHARTYPE>::ToString("")),
		minusActive_(true),
		isOk_(true),
		argv_from_cmdline_(nullptr),
		helpRequested_(false)
	{
		typename TSTR<CHARTYPE>::type s(c);
		// Build argc, argv from s. We must add a dummy first element for
		// progName because parse() expects it!!
		std::deque<typename TSTR<CHARTYPE>::type> args;
		//args.push_back(TSTR<CHARTYPE>::ToString(""));
		typename TSTRSTREAM<CHARTYPE>::I is(s);
		while (is.good())
		{
			typename TSTR<CHARTYPE>::type t;
			is>>t;
			args.push_back(t);
		}

		argv_from_cmdline_.reset(new PCHARTYPE[args.size()]);
		PCHARTYPE* p = argv_from_cmdline_.get();
		for (typename std::deque<typename TSTR<CHARTYPE>::type>::const_iterator
			iter = args.begin();
			iter != args.end();++iter)
		{
			*p++ = const_cast<CHARTYPE*>(iter->c_str());
		}

		//PCHARTYPE* pc = argv_from_cmdline_.get();
		parse(args.size(), argv_from_cmdline_.get());
	}

	template<typename CHARTYPE>
	inline void argstream<CHARTYPE>::parse(int argc, CHARTYPE const* const argv[])
	{
		// Run thru all arguments.
		// * it has -- in front : it is a long name option, if remainder is empty,
		//                        it is an error
		// * it has - in front  : it is a sequence of short name options, if
		//                        remainder is empty, deactivates option (- will
		//                        now be considered a _TCHAR).
		// * if any other _TCHAR, or if option was deactivated
		//                      : it is a value. Values are split in parameters
		//                      (immediately follow an option) and pure values.
		// Each time a value is parsed, if the previously parsed argument was an
		// option, then the option is linked to the value in case of it is a
		// option with parameter.  The subtle point is that when several options
		// are given with short names (ex: -abc equivalent to -a -b -c), the last
		// parsed option is -c).
		// Since we use map for option, any successive call overides the previous
		// one: foo -a -b -a hello is equivalent to foo -b -a hello
		// For values it is not true since we might have several times the same
		// value.
		value_iterator* lastOption = NULL;
		for (CHARTYPE** a = const_cast<CHARTYPE**>(argv),**astop=a+argc;++a!=astop;)
		{
			typename TSTR<CHARTYPE>::type s(*a);
			if (minusActive_ && s[0] == '-')
			{
				if (s.size() > 1 && s[1] == '-')
				{
					if (s.size() == 2)
					{
						minusActive_ = false;
						continue;
					}
					lastOption = &(options_[s.substr(2)] = values_.end());
				}
				else
				{
					if (s.size() > 1)
					{
						// Parse all _TCHARs, if it is a minus we have an error
						for (typename TSTR<CHARTYPE>::type::const_iterator cter = s.begin();
							++cter != s.end();)
						{
							if (*cter == '-')
							{
								isOk_ = false;
								typename TSTRSTREAM<CHARTYPE>::O os;
								os<<TSTR<CHARTYPE>::ToString("- in the middle of a switch ")<<a;
								errors_.push_back(os.str());
								break;
							}
							lastOption = &(options_[typename TSTR<CHARTYPE>::type(1,*cter)] = values_.end());
						}
					}
					else
					{
						isOk_ = false;
						errors_.push_back(TSTR<CHARTYPE>::ToString("Invalid argument -"));
						break;
					}
				}
			}
			else
			{
				values_.push_back(s);
				if (lastOption != NULL)
				{
					*lastOption = --values_.end();
				}
				lastOption = NULL;
			}
		}
#ifdef ARGSTREAM_DEBUG
		for (std::map<_TSTRING,value_iterator>::const_iterator
			iter = options_.begin();iter != options_.end();++iter)
		{
			TSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: option ") << iter->first;
			if (iter->second != values_.end())
			{
				TSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString(" -> ") << *(iter->second);
			}
			std::cout<<std::endl;
		}
		for (std::list<_TSTRING>::const_iterator
			iter = values_.begin();iter != values_.end();++iter)
		{
			TSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: value ") << *iter << std::endl;
		}
#endif // ARGSTREAM_DEBUG
	}

	template<typename CHARTYPE>
	inline bool
	argstream<CHARTYPE>::isOk() const
	{
		return isOk_;
	}

	template<typename CHARTYPE>
	inline bool
	argstream<CHARTYPE>::helpRequested() const
	{
		return helpRequested_;
	}

	template<typename CHARTYPE>
	inline typename TSTR<CHARTYPE>::type
	argstream<CHARTYPE>::usage() const
	{
		typename TSTRSTREAM<CHARTYPE>::O os;
		if (copyright_.size())
		{
			os << copyright_ << std::endl << std:: endl;
		}
		os<<TSTR<CHARTYPE>::ToString("Usage: ")<<progName_<<cmdLine_<<std::endl;
		unsigned int lmax = 0;
		for (typename std::deque<help_entry>::const_iterator iter = argHelps_.begin();
                     iter != argHelps_.end();++iter)
		{
			if (lmax<iter->first.size()) lmax = iter->first.size();
		}
		for (typename std::deque<help_entry>::const_iterator iter = argHelps_.begin();
                     iter != argHelps_.end();++iter)
		{
			os << '\t' << iter->first << typename TSTR<CHARTYPE>::type(lmax-iter->first.size(),' ')
				<< TSTR<CHARTYPE>::ToString(" : ") << iter->second << std::endl;
		}

		// Append the examples
		if (argExamples_.size())
		{
			int i = 1;
			for (auto example:argExamples_)
			{
				os << std::endl
					<< TSTR<CHARTYPE>::ToString("----------") << std::endl
					<< TSTR<CHARTYPE>::ToString("Example ") << i++ << std::endl
					<< TSTR<CHARTYPE>::ToString("----------") << std::endl
					<< TSTR<CHARTYPE>::ToString("Command: ") << std::endl
					<< TSTR<CHARTYPE>::ToString("\t") << example.first << std::endl << std::endl
					<< TSTR<CHARTYPE>::ToString("Description:") << std::endl
					<< TSTR<CHARTYPE>::ToString("\t") << example.second << std::endl;
			}
		}
		return os.str();
	}

	template<typename CHARTYPE>
	inline typename TSTR<CHARTYPE>::type
	argstream<CHARTYPE>::errorLog() const
	{
		typename TSTR<CHARTYPE>::type s;
		for(typename std::deque<typename TSTR<CHARTYPE>::type>::const_iterator iter = errors_.begin();
                    iter != errors_.end();++iter)
		{
			s += *iter;
			s += '\n';
		}
		return s;
	}

	template<typename CHARTYPE>
	inline RESULT_OF_PARSE
	argstream<CHARTYPE>::defaultErrorHandling(bool ignoreUnused) const
	{
		if (helpRequested_)
		{
			return RESULT_OF_PARSE::PARSED_ERR_HELP_REQUESTED;
		}
		if (!isOk_)
		{
			return RESULT_OF_PARSE::PARSED_ERR_OTHER;
		}
		if (!ignoreUnused &&
			(!values_.empty() || !options_.empty()))
		{
			return RESULT_OF_PARSE::PARSED_ERR_UNUSED_PARAMETER;
		}
		return RESULT_OF_PARSE::PARSED_OK;
	}

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Implementation of global functions
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	template <typename CHARTYPE, typename T>
	inline ValueHolder<CHARTYPE, T>
	parameter(
		CHARTYPE s,
		const CHARTYPE* l,
		T& b,
		const CHARTYPE* desc,
		bool mandatory)
	{
		return ValueHolder<CHARTYPE, T>(s,l,b,desc,mandatory);
	}

	/*
	template<typename CHARTYPE, typename T, typename O>
	inline ValuesHolder<CHARTYPE, T, O>
	values(
		const O& o,
		const CHARTYPE* desc,
		int len)
	{
		return ValuesHolder<CHARTYPE, T, O>(o, desc, len);
	}
	*/

	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>
	option(
		CHARTYPE s,
		const CHARTYPE* l,
		bool& b,
		const CHARTYPE* desc)
	{
		return OptionHolder<CHARTYPE>(s, l, b, desc);
	}
	template<typename CHARTYPE>
	inline OptionHolder<CHARTYPE>
	help()
	{
		assert(false); // Should not run the default path.
		return {};
	}
	template<>
	inline OptionHolder<wchar_t>
	help()
	{
		return OptionHolder<wchar_t>(L'h', L"Help", L"Display this help");
	}

	template<>
	inline OptionHolder<char>
	help()
	{
		return OptionHolder<char>('h', "Help", "Display this help");
	}

	template<typename CHARTYPE>
	inline ExampleHolder<CHARTYPE>
	example(const CHARTYPE* cmdline, const CHARTYPE* desc)
	{
		return ExampleHolder<CHARTYPE>(cmdline, desc);
	}

	template<typename CHARTYPE>
	inline CopyrightHolder<CHARTYPE>
	copyright(const CHARTYPE* copyright)
	{
		return CopyrightHolder<CHARTYPE>(copyright);
	}

	template<typename CHARTYPE, typename T>
	inline argstream<CHARTYPE>&
	operator >>(argstream<CHARTYPE>& s, ValueHolder<CHARTYPE, T> const& v)
	{
		// Search in the options if there is any such option defined either with a
		// short name or a long name. If both are found, only the last one is
		// used.
#ifdef ARGSTREAM_DEBUG
		TSTRSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: searching ")
			<< v.shortName_<<L" "<<v.longName_<<std::endl;
#endif
		typename argstream<CHARTYPE>::help_entry entry(v.name(), v.description());
		s.argHelps_.push_back(entry);
		if (v.mandatory_)
		{
			if (!v.shortName_.empty())
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" -");
				s.cmdLine_ += v.shortName_;
			}
			else
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" --");
				s.cmdLine_ += v.longName_;
			}
			s.cmdLine_ += TSTR<CHARTYPE>::ToString(" value");
		}
		else
		{
			if (!v.shortName_.empty())
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" [-");
				s.cmdLine_ += v.shortName_;
			}
			else
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" [--");
				s.cmdLine_ += v.longName_;
			}
			s.cmdLine_ += TSTR<CHARTYPE>::ToString(" value]");
		}
		typename std::map<typename TSTR<CHARTYPE>::type, typename argstream<CHARTYPE>::value_iterator>::iterator iter =
			s.options_.find(v.shortName_);
		if (iter == s.options_.end())
		{
			iter = s.options_.find(v.longName_);
		}
		if (iter != s.options_.end())
		{
			if (iter->second != s.values_.end())
			{
#ifdef ARGSTREAM_DEBUG
				TSTRSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: found value ")
					<< *(iter->second)<<std::endl;
#endif
				ValueParser<CHARTYPE, T> p;
				*(v.value_) = p(*(iter->second));
				// The option and its associated value are removed, the subtle thing
				// is that someother options might have this associated value too,
				// which we must invalidate.
				s.values_.erase(iter->second);

				/* Disabled by Levski Weng, if one item of the values_ is removed, you cannot compare them any more.
				for (std::map<typename TSTR<CHARTYPE>::type, typename argstream<CHARTYPE>::value_iterator>::iterator
					jter = s.options_.begin();jter != s.options_.end();++jter)
				{
					if (jter->second == iter->second)
					{
						jter->second = s.values_.end();
					}
				}
				*/
				s.options_.erase(iter);
			}
			else
			{
				s.isOk_ = false;
				typename TSTRSTREAM<CHARTYPE>::O os;
				os	<< TSTR<CHARTYPE>::ToString("No value following switch ") << iter->first
					<< TSTR<CHARTYPE>::ToString(" on command line");
				s.errors_.push_back(os.str());
			}
		}
		else
		{
			if (v.mandatory_)
			{
				s.isOk_ = false;
				typename TSTRSTREAM<CHARTYPE>::O os;
				os<< TSTR<CHARTYPE>::ToString("Mandatory parameter ");
				if (!v.shortName_.empty()) os<<'-'<<v.shortName_;
				if (!v.longName_.empty())
				{
					if (!v.shortName_.empty()) os<<'/';
					os << TSTR<CHARTYPE>::ToString("--") << v.longName_;
				}
				os<< TSTR<CHARTYPE>::ToString(" missing");
				s.errors_.push_back(os.str());
			}
		}
		return s;
	}

	template<typename CHARTYPE, typename T>
	inline argstream<CHARTYPE>&
	operator >>(argstream<CHARTYPE>& s, ValueHolder<CHARTYPE, bool> const& v)
	{
		// Search in the options if there is any such option defined either with a
		// short name or a long name. If both are found, only the last one is
		// used.
#ifdef ARGSTREAM_DEBUG
		TSTRSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: searching ")
			<< v.shortName_<<L" "<<v.longName_<<std::endl;
#endif
		s.argHelps_.push_back(argstream<CHARTYPE>::help_entry(v.name(),v.description()));
		if (v.mandatory_)
		{
			if (!v.shortName_.empty())
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" -");
				s.cmdLine_ += v.shortName_;
			}
			else
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" --");
				s.cmdLine_ += v.longName_;
			}
		}
		else
		{
			if (!v.shortName_.empty())
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" [-");
				s.cmdLine_ += v.shortName_;
			}
			else
			{
				s.cmdLine_ += TSTR<CHARTYPE>::ToString(" [--");
				s.cmdLine_ += v.longName_;
			}
			s.cmdLine_ += TSTR<CHARTYPE>::ToString("]");

		}
		typename std::map<typename TSTR<CHARTYPE>::type, typename argstream<CHARTYPE>::value_iterator>::iterator iter =
			s.options_.find(v.shortName_);
		if (iter == s.options_.end())
		{
			iter = s.options_.find(v.longName_);
		}
		if (iter != s.options_.end())
		{
			if (iter->second != s.values_.end())
			{
#ifdef ARGSTREAM_DEBUG
				TSTRSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: found value ")
					<< *(iter->second)<<std::endl;
#endif
				ValueParser<CHARTYPE, T> p;
				*(v.value_) = p(*(iter->second));
				// The option and its associated value are removed, the subtle thing
				// is that someother options might have this associated value too,
				// which we must invalidate.
				// Modified by Levski Weng
				//s.values_.erase(iter->second);
				for (typename std::map<typename TSTR<CHARTYPE>::type, typename argstream<CHARTYPE>::value_iterator>::iterator
					jter = s.options_.begin();jter != s.options_.end();++jter)
				{
					if (jter->second == iter->second)
					{
						jter->second = s.values_.end();
					}
				}
				s.options_.erase(iter);
			}
			else
			{
				*(v.value_) = true;
			}
		}
		else
		{
			if (v.mandatory_)
			{
				s.isOk_ = false;
				typename TSTRSTREAM<CHARTYPE>::O os;
				os<< TSTR<CHARTYPE>::ToString("Mandatory parameter ");
				if (!v.shortName_.empty()) os<<'-'<<v.shortName_;
				if (!v.longName_.empty())
				{
					if (!v.shortName_.empty()) os<<'/';
					os << TSTR<CHARTYPE>::ToString("--") << v.longName_;
				}
				os<< TSTR<CHARTYPE>::ToString(" missing");
				s.errors_.push_back(os.str());
			}
		}
		return s;
	}


	template<typename CHARTYPE>
	inline argstream<CHARTYPE>&
	operator >>(argstream<CHARTYPE>& s, OptionHolder<CHARTYPE> const& v)
	{
		// Search in the options if there is any such option defined either with a
		// short name or a long name. If both are found, only the last one is
		// used.
#ifdef ARGSTREAM_DEBUG
		TSTRSTREAM<CHARTYPE>::COUT << TSTR<CHARTYPE>::ToString("DEBUG: found value ")
			<< v.shortName_<<L" "<<v.longName_<<std::endl;
#endif
		typename argstream<CHARTYPE>::help_entry entry(v.name(), v.description());
		s.argHelps_.push_back(entry);
		{
			typename TSTR<CHARTYPE>::type c;
			if (!v.shortName_.empty())
			{
				c += TSTR<CHARTYPE>::ToString(" [-");
				c += v.shortName_;
			}
			else
			{
				c += TSTR<CHARTYPE>::ToString(" [--");
				c += v.longName_;
			}
			c += TSTR<CHARTYPE>::ToString("]");
			s.cmdLine_ = c+s.cmdLine_;
		}

		if (s.options_.find(TSTR<CHARTYPE>::ToString('h')) != s.options_.end() ||
			s.options_.find(TSTR<CHARTYPE>::ToString("help")) != s.options_.end() )
		{
			s.helpRequested_ = true;
		}
		typename std::map<typename TSTR<CHARTYPE>::type, typename argstream<CHARTYPE>::value_iterator>::iterator iter =
			s.options_.find(v.shortName_);
		if (iter == s.options_.end())
		{
			iter = s.options_.find(v.longName_);
		}
		if (iter != s.options_.end())
		{
			// If we find counterpart for value holder on command line then the
			// option is true and if an associated value was found, it is ignored
			if (v.value_ != NULL)
			{
				*(v.value_) = true;
			}
			// The option only is removed
			s.options_.erase(iter);
		}
		else
		{
			if (v.value_ != NULL)
			{
				*(v.value_) = false;
			}
		}
		return s;
	}

};
#endif // ARGSTREAM_H