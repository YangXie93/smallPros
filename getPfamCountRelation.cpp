#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <list>
#include<iterator>

using namespace Rcpp;

//[[Rcpp::export]]

List getPfamCountRelation(List refList, std::vector<int> orgVec, std::vector<int> query){
    
    std::vector<int> missing (orgVec.size(),0);
    std::vector<int> exceeding (orgVec.size(),0);
    std::vector<int> total (orgVec.size(),0);
    std::vector<bool> isIn (orgVec.size(),true);
    
    std::vector<int>::iterator missingIt = missing.begin();
    std::vector<int>::iterator exceedingIt = exceeding.begin();
    std::vector<int>::iterator totalIt = total.begin();
    std::vector<bool>::iterator isInIt = isIn.begin();
    
    std::vector<int>::iterator queryIt = query.begin();
    
    List::iterator orgListIt;
    std::vector<int>::iterator pfamCountVecIt;
    std::vector<int>::iterator orgListVecIt;
    
    int i = 0;
    
    for(List::iterator refListIt = refList.begin(); refListIt != refList.end(); refListIt++){
        
        if(*(refListIt) == NULL){
            if(*next(queryIt,i) != 0){
                for(int j =0 ;j < exceeding.size();j++){
                    *next(exceedingIt,j) += *next(queryIt,i);
                }
            }
        }
        else{
            pfamCountVecIt = (*refListIt).begin();
            orgListIt = (*next(refListIt.begin())).begin();
            int tmp;
            
            for(int j = 0; j < (*refListIt).size();j++){
                
                for(int n = 0; n < (*orgListIt).size();n++){
                    
                    tmp = (*pfamCountVecIt)[n] - *next(queryIt,i);
                    
                    *next(totalIt,(*orgListIt)[n]) += (*pfamCountVecIt)[n];
                    if(tmp < 0){
                        *next(exceedingIt,(*orgListIt)[n]) += -1* tmp;
                    }
                    else{
                        *next(missingIt,(*orgListIt)[n]) += tmp;
                    }
                    
                }
                
                
                pfamCountVecIt++;
                orgListIt++;
            }
            
        }
       
       
       i++; 
    }
    
    List res = List::create(Named("total") = total,Named("missing") = missing,Named("exceeding") = exceeding);
    
    return res;
}
