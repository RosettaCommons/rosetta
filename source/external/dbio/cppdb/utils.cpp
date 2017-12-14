///////////////////////////////////////////////////////////////////////////////
//                                                                             
//  Copyright (C) 2010-2011  Artyom Beilis (Tonkikh) <artyomtnk@yahoo.com>     
//                                                                             
//  Distributed under:
//
//                   the Boost Software License, Version 1.0.
//              (See accompanying file LICENSE_1_0.txt or copy at 
//                     http://www.boost.org/LICENSE_1_0.txt)
//
//  or (at your opinion) under:
//
//                               The MIT License
//                 (See accompanying file MIT.txt or a copy at
//              http://www.opensource.org/licenses/mit-license.php)
//
///////////////////////////////////////////////////////////////////////////////
#define CPPDB_SOURCE
#include <cppdb/utils.h>
#include <cppdb/errors.h>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <locale>

namespace cppdb {
	std::string format_time(std::tm const &v)
	{
		char buf[64]= {0};
		strftime(buf,sizeof(buf),"%Y-%m-%d %H:%M:%S",&v);
		return std::string(buf);
	}

	std::tm parse_time(std::string const &v)
	{
		if(strlen(v.c_str())!=v.size())
			throw bad_value_cast();
		return parse_time(v.c_str());
	}
	std::tm parse_time(char const *v)
	{
		std::tm t=std::tm();
		int n;
		double sec = 0;
		n = sscanf(v,"%d-%d-%d %d:%d:%lf",
			&t.tm_year,&t.tm_mon,&t.tm_mday,
			&t.tm_hour,&t.tm_min,&sec);
		if(n!=3 && n!=6) 
		{
			throw bad_value_cast();
		}
		t.tm_year-=1900;
		t.tm_mon-=1;
		t.tm_isdst = -1;
		t.tm_sec=static_cast<int>(sec);
		if(mktime(&t)==-1)
			throw bad_value_cast();
		return t;
	}

	namespace {
		bool is_blank_char(char c)
		{
			return c==' ' || c=='\t' || c=='\r' || c=='\n' || c=='\f';
		}
		std::string trim(std::string const &s)
		{
			if(s.empty())
				return s;
			size_t start=0,end=s.size()-1;
			while(start < s.size() && is_blank_char(s[start])) {
				start++;
			}
			while(end > start && is_blank_char(s[end])) {
				end--;
			}
			return s.substr(start,end-start+1);
		}
	}

	void parse_connection_string(	std::string const &connection_string,
					std::string &driver,
					std::map<std::string,std::string> &params)
	{
		params.clear();
		size_t p = connection_string.find(':');
		if( p == std::string::npos )
			throw cppdb_error("cppdb::Invalid connection string - no driver given");
		driver = connection_string.substr(0,p);
		p++;
		while(p<connection_string.size()) {
			size_t n=connection_string.find('=',p);
			if(n==std::string::npos)
				throw cppdb_error("Invalid connection string - invalid property");
			std::string key = trim(connection_string.substr(p,n-p));
			p=n+1;
			std::string value;
			while(p<connection_string.size() && is_blank_char(connection_string[p]))
			{
				++p;
			}
			if(p>=connection_string.size()) {
				/// Nothing - empty property
			}
			else if(connection_string[p]=='\'') {
				p++;
				while(true) {
					if(p>=connection_string.size()) {
						throw cppdb_error("Invalid connection string unterminated string");
					}
					if(connection_string[p]=='\'') {
						if(p+1 < connection_string.size() && connection_string[p+1]=='\'') {
							value+='\'';
							p+=2;
						}
						else {
							p++;
							break;
						}
					}
					else {
						value+=connection_string[p];
						p++;
					}
				}
			}
			else {
				size_t n=connection_string.find(';',p);
				if(n==std::string::npos) {
					value=trim(connection_string.substr(p));
					p=connection_string.size();
				}
				else {
					value=trim(connection_string.substr(p,n-p));
					p=n;
				}
			}
			if(params.find(key)!=params.end()) {
				throw cppdb_error("cppdb::invalid connection string duplicate key");
			}
			params[key]=value;
			while(p<connection_string.size()) {
				char c=connection_string[p];
				if(is_blank_char(c))
					++p;
				else if(c==';') {
					++p;
					break;
				}
			}
		}
	} //
	std::string connection_info::get(std::string const &prop,std::string const &default_value) const
	{
		auto p=properties.find(prop);
		if(p==properties.end())
			return default_value;
		else
			return p->second;
	}

	bool connection_info::has(std::string const &prop) const
	{
		return properties.find(prop) != properties.end();
	}

	int connection_info::get(std::string const &prop,int default_value) const
	{
		auto p=properties.find(prop);
		if(p==properties.end())
			return default_value;
		std::istringstream ss;
		ss.imbue(std::locale::classic());
		ss.str(p->second);
		int val;
		ss >> val;
		if(!ss || !ss.eof()) {
			throw cppdb_error("cppdb::connection_info property " + prop + " expected to be integer value");
		}
		return val;
	}

}
