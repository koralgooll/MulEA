#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world_1() {
	CharacterVector x = CharacterVector::create("foo", "bar");
	NumericVector y = NumericVector::create( 0.0, 1.0 ) ;
	List z = List::create( x, y ) ;
	return z ;
}

// [[Rcpp::export]]
List rcpp_hello_world_2() {
  CharacterVector x = CharacterVector::create("twoo", "twoooo");
  NumericVector y = NumericVector::create( 2.0, 2.0 ) ;
  List z = List::create( x, y ) ;
  return z ;
}
