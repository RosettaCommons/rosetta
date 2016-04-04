#ifndef self_test_01_hpp
#define self_test_01_hpp

// #include <string>
#include <vector>
#include <map>
#include <iostream>
//#include <memory>
// #include <complex>
//#include <utility>
//#include <iostream>

#include <self_test.incl.hpp>

namespace aa {

// class AA {
// public:
// 	int n;
// 	float f;
// 	double d;
// };

// typedef int aa_INT;

// inline void foo_aa(A& a, aa_INT) {}

// typedef void (* UtilityExitCallBack)(void);

// void set_main_exit_callback( UtilityExitCallBack = 0 );

}

//std::string foo1(std::string s = std::string("aaa") ) { return s; }

//std::string foo2(std::string s = "aaa" ) { return s; }
//void foo3(A a=A() ) {}

// inline bool foo(bool a, bool&b, aa::A) { return true; }
// std::string foo(int a) { return "foo1"; }
// std::string foo(int a, int b) { return "foo2"; }
// std::string foo(int a, int b, int c) { return "foo3: a="+std::to_string(a)+" b="+std::to_string(b)+" c="+std::to_string(c); }


// inline void foo_A(aa::A& a, aa::aa_INT) {}

enum E1 {E1A, E1B};


namespace utility {

// void foo_p(std::pair<int, int> &) {}
// void foo_cp(std::pair<int, int> const &) {}
// void foo_t(std::tuple<int, int> &) {}
// void foo_ct(std::tuple<int, int> const &) {}
//void foo_a(aa::AA const &) {}
//void foo_i(int const &) {}

// // template <typename T>
// class A {
// };
// // void foo_t(A<enum aaaa::EA> &) {}

// struct no_property {};
// template <class Tag, class T, class Base = no_property>
// struct property {
// };

// template <class T> void foo_t(T &) {}


// template <> void foo_t(property<aaaa::C::QWE, int> &) {};
//void foo_t(property<enum aaaa::EA, int> &) {}



class A {};

void foo(int a = 1) {}


// std::vector<int> foo_v(std::vector<std::string> &)
// {
// 	std::vector<int> v {1,2,3};
// 	return v;
// }

// std::vector<A> foo_a()
// {
// 	std::vector<A> v {A(),A(),A()};
// 	return v;
// }

// void foo_at(std::vector<A>) {}

// std::map<int, int> foo_m()
// {
// 	std::map<int, int> m {{1,1}, {2,3} };
// 	return m;
// }

// class A //: public std::enable_shared_from_this<A>
// {
// 	~A() {}
// public:
// 	A(int &r) : c(r) {}

// 	static void delete_(A* p) { p->~A(); }

// 	int a, b;
// 	int &c;
// };

// Test for binding data member with deleted operator=

// class MyData
// {
// public:
// 	MyData& operator=(const MyData& other) = delete;
// };

// typedef void * VoidP;

// class A
// {
// 	struct PR {};
// public:
// 	int a;
// 	// MyData m;

// 	void foo_pr(PR);
// 	// void foo_arg(aa::AA) {}
// 	// void foo_map(std::map<int, int>) {}

// 	//void *foo_m() { return nullptr; };
// 	void foo_m(VoidP) {};
// };

// template <typename T>
// class TT {
// public:
// 	explicit TT(int ) {}
// 	//void foo() {}
// };

// void foo_i(int) {}
// void foo_c(char) {}


// void foo_a(TT<aa::AA>) {}
// void foo_arg(aa::AA) {}

// template class TT<const aa::AA &>;
// template class TT<aa::AA &>;
// template class TT<aa::AA const>;

// template class TT<A::PR &>;

// void foo(std::_Ios_Openmode) {}

// struct VoidPointer
// {
// 	void *p;
// 	//explicit operator void*() const { return p; }
// };
// void * foo_null() { return nullptr; }
// void * foo_null1() { return (void*)1; }
// void foo_null2(void *p) { std::cout << p << std::endl; }



// template< typename T >
// class AtomID_Map
// {
// public:
// 	explicit AtomID_Map(T const & default_value_a ) {}

// 	//explicit AtomID_Map(AtomID_Map const & ) = default;
// 	//AtomID_Map & operator= (AtomID_Map const & ) = default;
// };

// // class R : public std::enable_shared_from_this<R>
// // {
// // public:
// // 	void foo_t(AtomID_Map<bool>) {}
// // };
// void foo_t(AtomID_Map<bool>) {}

//void foo_p(std::pair<int, int> ) {}

//void foo_complex(std::complex<double> a) {}

// struct S
// {
// 	S() {}
// 	S(aaaa::A &){}

// 	int a;
// 	double b;
// 	int c[10];
// };

// template<
// 	template< typename > class Array
// >
// class FArrayInitializer
// {};

// void foo_e( FArrayInitializer< aaaa::AT > ) {}


// class DP
// {
// public:
// 	DP(const double *){}
// };


// void foo_cdp(double *) {}
// void foo_cdpc(const double *) {}
// void foo_cd(const double) {}
// void foo_d(double) {}

//void foo_skip() {};

// template <class T> class N
// {
// 	enum OptionTypes {
// 		UNKNOWN_OPTION,
// 	};

// 	void private_foo(std::vector<OptionTypes> );

// public:
// 	void foo() {}
// };

// void foo_N(N<int>) {}


//void foo() {}
//void foo_string_const(std::string) {}

// void foo_string_const(std::string const &) {}
// void foo_string(std::string  &) {}
// void foo_int_const(int const &) {}
// void foo_int(int &) {}

// class CA {};
// template<typename T>
// void foo_aaaa(int) {}

// //template<>
// //void foo_aaaa<aaaa::A>(int);

// void foo()
// {
// 	foo_aaaa<aaaa::A>(0);
// }



// template<typename T>
// class TA
// {
// public:
// 	T value;
// };


// template<int S>
// class TI
// {
// public:
// 	int values[S];
// };


// class B : public TA<aaaa::A>
// {
// public:
// 	enum Color {Red, Blue};

// 	void foo(B::Color c = B::Color::Blue) {}

// 	int a;

// 	void foo_open_mode(std::ios_base::openmode open_mode = std::ios_base::out) {}

// 	/// @brief Open a file
// 	void
// 	open(
// 		std::string const & filename_a,
// 		std::ios_base::openmode open_mode = std::ios_base::out
// 	);

// 	void foo_m(B&&) {}
// };

// void foo_c(B::Color c = utility::B::Blue) {}

//void foo_A( std::shared_ptr< const aaaa::A > &) {}
//void foo_A( std::shared_ptr< TA<aaaa::A> > &) {}

// typedef TA<aaaa::A> TA_A;

// void foo_TA(TA<aaaa::A>) {}

// template<int T>
// void foo_T(TI<T> &) {}


// namespace sql_database {

// struct TransactionMode {
// 	enum e {
// 		none = 1,
// 		standard,
// 		chunk
// 	};
// };

// struct DatabaseMode {
// 	enum e {
// 		sqlite3 = 1,
// 		mysql,
// 		postgres
// 	};
// };

// TransactionMode::e
// transaction_mode_from_name(
// 	std::string transaction_mode);

// std::string
// name_from_transaction_mode(
// 	TransactionMode::e transaction_mode);

// DatabaseMode::e
// database_mode_from_name(
// 	std::string database_mode);

// std::string
// name_from_database_mode(
// 	DatabaseMode::e database_mode);

// }


//void foo_c(B::Color c = B::Red) {}

// typedef TA<aaaa::A> TA_A;

// class TC : public TA<aaaa::A> {};

// class A {
// public:
// 	int n;

// 	void foo() {};
// };

// template <typename T>
// class TUA;

// template <typename T>
// class TUA {
// public:
// 	T t;
// };


// typedef TUA<int> TUA_int;


// void foo_TUA(TUA_int a) {};


// template< bool >
// struct vectorL_IndexSelector
// {
// };


// /// @brief vectorL index type selector: Negative lower index specialization
// template<>
// struct vectorL_IndexSelector< false >
// {
// 	inline
// 	static
// 	bool
// 	ge() { return true; }
// };


// template <class T>
// void swap(TUA<T>) {}

// void foo_i(int **a) {}
// std::string foo_s(std::string const &a="aaa") { return a; }


// class TUAC : public TUA_int
// {};

// void foo(TUA_int) {}
// void foo(TUA<float>) {}

// namespace csi {

// enum E2 {E2A, E2B};

// template<class T>
// void foo_T() {}

// template void foo_T<int>();
// template void foo_T<char>();

// namespace inner1 {
// namespace inner2 {
// namespace inner3 {
// void foo() {}
// }
// }
// }

// class CSI_Sequence
// {
// public:
// 	/// @brief constructor
// 	CSI_Sequence(std::string sequence) {}

// 	operator std::string() const { return sequence_; }

// 	/// @brief operator to output our sequence so we can write: std::cout << CSI_SequenceObject
// 	//friend std::ostream & operator << (std::ostream & os, CSI_Sequence const &sq) { os << sq.sequence_; return os; }

// 	void foo(int) const {};
// private:
// 	std::string sequence_;
// };

// void foo_csi(char * const argv[]) {}

// class ESFT : public std::enable_shared_from_this<ESFT>
// {
// public:
// 	int a;
// };


// class ESFT2 : public ESFT
// {
// public:
// 	std::shared_ptr<ESFT2> shared_from_this() {
//         return shared_from_this();
//     }

// 	void * fn;
// 	int a;

// };

// void foo(std::string a = "abc") {}
// std:: string foo(int a) { return ""; }

//} // namespace csi

} // namespace utility

#endif // self_test_01_hpp
