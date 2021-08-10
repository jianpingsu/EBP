#include"BiljectivePara\BiljectivePara.h"
#include<string>
#include<iostream>
#include<fstream>
#include<io.h>
using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cerr << "Syntax: " << argv[0] << " <input mesh>" << endl;
		return -1;
	}
	if (2 == argc)
	{
		const string input_mesh = argv[1];
		
		BiljectivePara * bil_para = new BiljectivePara(input_mesh);
		
		bil_para->load();
		bil_para->parameterization();
		
		delete bil_para;
		bil_para = NULL;

		system("pause");
		return 0;
	}
}
