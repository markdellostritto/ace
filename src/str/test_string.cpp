// c++
#include <iostream>
// str
#include "str/print.hpp"
#include "str/string.hpp"

void test_number(){
    const int nstr=7;
    std::vector<std::string> svec(nstr);
    svec[0]="1234";
    svec[1]="-1234";
    svec[2]="1.234";
    svec[3]="1.234e10";
    svec[4]="1.234E10";
    svec[5]="1.234a10";
    svec[6]="Trust me, I'm a number.";
    //print
    char* str=new char[print::len_buf];
    std::cout<<print::buf(str,print::char_buf)<<"\n";
    std::cout<<print::title("TEST - NUMBER",str)<<"\n";
    for(int i=0; i<nstr; ++i){
        std::cout<<"num "<<string::number(svec[i].c_str())<<" str "<<svec[i]<<"\n";
    }
    std::cout<<print::buf(str,print::char_buf)<<"\n";
    delete[] str;
}

int main(int argc, char* argv[]){
	
    test_number();

	return 0;
}