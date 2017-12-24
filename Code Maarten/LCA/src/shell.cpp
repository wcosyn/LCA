#include "shell.h"

#include <vector>
using std::vector;

int Shell::shellsinit = 0;
vector< Shell* > Shell::shellsN= vector< Shell* >(22);
vector< Shell* > Shell::shellsP= vector< Shell* >(22);

vector< Shell > Shell::shells={
    Shell(0,0,1), //2

    Shell(0,1,3), //6
    Shell(0,1,1), //8
    
    Shell(0,2,5), //14
    Shell(1,0,1), //16
    Shell(0,2,3), //20

    Shell(0,3,7), //28

    Shell(1,1,3), //32
    Shell(0,3,5), //38
    Shell(1,1,1), //40
    Shell(0,4,9), //50

    Shell(0,4,7), //58
    Shell(1,2,5), //64
    Shell(1,2,3), //68
    Shell(2,0,1), //70
    Shell(0,5,11), //82

    Shell(0,5,9), //92
    Shell(1,3,7), //100
    Shell(1,3,5), //106
    Shell(2,1,3), //110
    Shell(2,1,1), //112
    Shell(0,6,13), //126
 };

void Shell::initializeShells()
{
    if(!shellsinit ) {
        shellsP[0]  = new Shell(0,0,1); //2

        shellsP[1]  = new Shell(0,1,3); //6
        shellsP[2]  = new Shell(0,1,1); //8

        shellsP[3]  = new Shell(0,2,5); //14
        shellsP[4]  = new Shell(1,0,1); //16
        shellsP[5]  = new Shell(0,2,3); //20

        shellsP[6]  = new Shell(0,3,7); //28

        shellsP[7]  = new Shell(1,1,3); //32
        shellsP[8]  = new Shell(0,3,5); //38
        shellsP[9]  = new Shell(1,1,1); //40
        shellsP[10] = new Shell(0,4,9); //50

        shellsP[11] = new Shell(0,4,7); //58
        shellsP[12] = new Shell(1,2,5); //64
        shellsP[13] = new Shell(1,2,3); //68
        shellsP[14] = new Shell(2,0,1); //70
        shellsP[15] = new Shell(0,5,11); //82

        shellsP[16] = new Shell(0,5,9); //92
        shellsP[17] = new Shell(1,3,7); //100
        shellsP[18] = new Shell(1,3,5); //106
        shellsP[19] = new Shell(2,1,3); //110
        shellsP[20] = new Shell(2,1,1); //112
        shellsP[21] = new Shell(0,6,13); //126

        shellsN[0]  = new Shell(0,0,1); //2
        shellsN[1]  = new Shell(0,1,3); //6
        shellsN[2]  = new Shell(0,1,1); //8
        shellsN[3]  = new Shell(0,2,5); //14
        shellsN[4]  = new Shell(1,0,1); //16
        shellsN[5]  = new Shell(0,2,3); //20

        shellsN[6]  = new Shell(0,3,7); //28

        shellsN[7]  = new Shell(1,1,3); //32
        shellsN[8]  = new Shell(0,3,5); //38
        shellsN[9]  = new Shell(1,1,1); //40
        shellsN[10] = new Shell(0,4,9); //50

        shellsN[11] = new Shell(0,4,7); //58
        shellsN[12] = new Shell(1,2,5); //64
        shellsN[13] = new Shell(1,2,3); //78
        shellsN[14] = new Shell(2,0,1); //70
        shellsN[15] = new Shell(0,5,11);//82

        shellsN[16] = new Shell(0,5,9); //92
        shellsN[17] = new Shell(1,3,7); //100
        shellsN[18] = new Shell(1,3,5); //106
        shellsN[19] = new Shell(2,1,3); //110
        shellsN[20] = new Shell(2,1,1); //112

        shellsN[21] = new Shell(0,6,13); //126

        shellsinit++;
    } else {
        shellsinit++;
    }
}

void Shell::deleteShells()
{
    shellsinit--;
    if( !shellsinit ) {
        for( unsigned int i= 0; i< shellsN.size(); i++) {
            delete shellsN[i];
        }
        for( unsigned int i= 0; i< shellsP.size(); i++) {
            delete shellsP[i];
        }
    }
}

// Also see shell.h
void Shell::get_shell_max( const int A,  int& shell_max,  int& max )
{
    if( A <= 2 ) {
        shell_max= 0;
        max= 2;
    } else if( A <= 6 ) {
        shell_max= 1;
        max= 6;
    } else if( A <= 8 ) {
        shell_max= 2;
        max= 8;
    } else if( A <= 14 ) {
        shell_max= 3;
        max= 14;
    } else if( A <= 16 ) {
        shell_max= 4;
        max= 16;
    } else if( A <= 20) {
        shell_max= 5;
        max= 20;
    } else if( A <= 28 ) {
        shell_max= 6;
        max= 28;
    } else if( A <= 32 ) {
        shell_max= 7;
        max= 32;
    } else if( A <= 38 ) {
        shell_max= 8;
        max= 38;
    } else if( A <= 40 ) {
        shell_max= 9;
        max= 40;
    } else if( A <= 50 ) {
        shell_max= 10;
        max= 50;
    } else if( A <= 58 ) {
        shell_max= 11;
        max= 58;
    } else if( A <= 64 ) {
        shell_max= 12;
        max= 64;
    } else if( A <= 68 ) {
        shell_max= 13;
        max= 68;
    } else if( A <= 70 ) {
        shell_max= 14;
        max= 70;
    } else if( A <= 82 ) {
        shell_max= 15;
        max= 82;
    } else if( A <= 92 ) {
        shell_max= 16;
        max= 92;
    } else if( A <= 100 ) {
        shell_max= 17;
        max= 100;
    } else if( A <= 106 ) {
        shell_max= 18;
        max= 106;
    } else if( A <= 110 ) {
        shell_max= 19;
        max= 110;
    } else if( A <= 112 ) {
        shell_max= 20;
        max= 112;
    } else if( A <= 126 ) {
        shell_max= 21;
        max= 126;
    }
}

