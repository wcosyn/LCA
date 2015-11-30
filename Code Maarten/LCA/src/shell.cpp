#include "shell.h"

int Shell::shellsinit = 0;
vector< Shell* > Shell::shellsN= vector< Shell* >(22);
vector< Shell* > Shell::shellsP= vector< Shell* >(22);

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
