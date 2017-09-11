library(Rcpp) 
# required for inline Rcpp calls 
library(inline)   
# write the C++ code 
do_rcpp_src <- ' 
// get data from the input data frame 
Rcpp::DataFrame cohort(the_cohort);\n 
// now extract columns by name from 
// the data fame into C++ vectors 
std::vector<double> age_v = Rcpp::as< std::vector<double> >(cohort["age"]);
std::vector<int> female_v = Rcpp::as< std::vector<int> >(cohort["female"]); 
std::vector<int> ily_v = Rcpp::as< std::vector<int> >(cohort["ily"]);
// create a new variable v_prob for export 
std::vector<double> v_prob (ily_v.size());
// iterate over data frame to calculate v_prob 
for (int i = 0; i < v_prob.size() ;  i++) { 
  v_prob[i] = vaccinate_cxx(age_v[i],female_v[i],ily_v[i]); 
}   
// export the old with the new in a combined data frame 
return 
Rcpp::DataFrame::create( Named("age")= age_v, Named("female") = female_v, Named("ily") = ily_v, Named("p") = v_prob);
'
# write the helper function also in C++ 
vaccinate_cxx_src <- ' 
using namespace std;
double vaccinate_cxx (double age, int female, int ily){ 
// this is based on some pretend regression equation 
double p = 0.25 + 0.3 * 1/(1-exp(0.004*age)) + 0.1 *ily;
// use some logic 
p = p * (female ? 1.25 : 0.75);
// boundary checking 
p = max(0.0,p);
p = min(1.0,p);
return(p);
} ' 
# create an R function to call the C++ code 
do_rcpp <- cxxfunction(signature(the_cohort="data.frame"), 
                       do_rcpp_src, plugin="Rcpp", 
                       includes=c('#include <cmath>', vaccinate_cxx_src)) 
#- See more at: http:\n//www.babelgraph.org/wp/?p=358#sthash.8bIJrV1p.dpuf

############


