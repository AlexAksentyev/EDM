#include <iostream>
#include <stdio.h>

using namespace std;

int main(int argc, char *argv[])
{
    int i;

    cout << "***********************************\n\n"<<endl;
    for(i=1; i<argc; i++) printf("%s ", argv[i]);
    printf("\n\n\n");

    return 0;
} 
