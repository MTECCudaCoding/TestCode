/***********************************************************************
 ***********************************************************************
 **********************Coincidencias Similares**************************
 ***********************************************************************
 ****************Analizador de un día contra un día*********************
 ***********************************************************************
 ***********************************************************************
 ***********************************************************************
 ****************************Diego GR***********************************
 ***********************************************************************
*/
/*
cd /home/monitec/Documentos/ProgramasEnC/ExtractorSimilares/
gcc -o a AnalisadorMuestras.c
./a '03-Costa Rica Nacional-20151031000000-50-0000-canal 7-0000-Servidor15-0000-0000.DAT'
./a '01-Costarica05-20151022071837-600-90.3 MHz-7E4885960-100-PCx_10_14-2-MNTC0023.dat'
 *
 * */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#define CANTIDAD_MAXIMA_RECONOCIMIENTOS 700
#define NUM_BYTES_FRECUENCY 10
#define NUM_FREC_COM 4
#define LEN_FREC_COM 32
#define NUM_INT_ESCALON 4
#define TAM_ESCALON 128



#define PORCENTAJE_BARRA_MINIMO_PERMITIDO 0.35f
#define PORCENTAJE_MINIMO_COINCIDENCIAS_PERMITIDO 0.8f

#define PUNTOS_ANALISIS_PRIMER_FILTRO 32
#define TAMANO_MINIMO_MUESTRA_SEGUNDOS 10
#define RESTRICCION_CERCANIA_SEGUNDOS 10

#define NUM_BYTES_SHORT 2
const float CONST_MEAN[] = {0.973f,0.93f,0.905f,0.9f,0};

typedef struct{
    unsigned int nums[NUM_INT_ESCALON];
    char bloqueado;
    char vacio;
}entero128;

typedef struct{
    unsigned int inicio;
    unsigned int duracion;
}reconocimiento;
typedef struct{
    unsigned int indice;
    unsigned int duracion;
    float porcentaje;
}coincidencia;

typedef struct{
	unsigned int inicio;
	unsigned int tamano;
	int cantidadCoincidencias;
	int cantidadReconocimientos;
	coincidencia coincidencias[CANTIDAD_MAXIMA_RECONOCIMIENTOS];
	reconocimiento reconocimientos[CANTIDAD_MAXIMA_RECONOCIMIENTOS];
}huella;


unsigned int multiplicadorTiempo = 1;
unsigned int SEGUNDOS_POR_DIA  = 86400;
void analizar(const char * pRutaDat);
int fsize(const char *filename);
void leerSuperDat(const char * pPathDatIn, int lenBuffer,unsigned short **  sliceDAT);
void discretizarDAT(const int pLenBuffer, unsigned short **pSliceDAT, entero128 *pMatrizOUT);
void extraerMuestra(int pIndiceInicio, const int pTamano, entero128 *pMatrizOUT_DAT,  huella *pHuella_Muestra);
void revisarHuella(int pIndiceInicio, const int pLenMatriz, entero128 *pMatrizOUT_DAT, huella * pHuella_Muestra, const int pLimite, const int pAvance, const float pPORCENTAJE, unsigned int * pArregloHamming);
void obtenerHamming(int pIndiceInicio, const int pLenMatriz, entero128 *pMatrizOUT_DAT, int pTam_Muestra, int pIni_Muestra, unsigned int * pArregloHamming,const int pAvance);


int main( int arc, char **argv ){

	analizar(argv[1]);
	return 0;
}

void analizar(const char * pRutaDat){
	unsigned int * GLOBAL_ARREGLO_HAMMING_h;
	entero128 * GLOBAL_MATRIZ_DAT_h;
	int lenBuffer = fsize(pRutaDat)/NUM_BYTES_FRECUENCY;


	if(lenBuffer <= 126000){
		SEGUNDOS_POR_DIA = 630;
	}
	int vEscalonAnalisis = 0, vIndMatriz, vTamMuestraDAT, vTamMuestraMatriz, vCantidadHuellas = 0, vIndCoincidencias = 0;
	multiplicadorTiempo = lenBuffer/SEGUNDOS_POR_DIA;

    unsigned short **sliceDAT;
    sliceDAT = (unsigned short **)malloc(NUM_FREC_COM*sizeof(short*));
    int vIndFrecuency;
    for(vIndFrecuency = 0;vIndFrecuency<NUM_FREC_COM;vIndFrecuency++){
        sliceDAT[vIndFrecuency] = (unsigned short *)malloc(lenBuffer*sizeof(short));
	}
	leerSuperDat(pRutaDat, lenBuffer, sliceDAT);

	int vlenMatriz = lenBuffer - LEN_FREC_COM;

	GLOBAL_ARREGLO_HAMMING_h = (unsigned int *) malloc(vlenMatriz*sizeof(unsigned int));

	GLOBAL_MATRIZ_DAT_h = (entero128 *)malloc(vlenMatriz*sizeof(entero128));
/*
	for(vIndMatriz = 0; vIndMatriz < vlenMatriz; vIndMatriz ++){
		for(vEscalonAnalisis = 0; vEscalonAnalisis < NUM_INT_ESCALON; vEscalonAnalisis ++){
			GLOBAL_MATRIZ_DAT_h[vIndMatriz].nums[vEscalonAnalisis] = 15000;
		}
	}
*/

	discretizarDAT(lenBuffer, sliceDAT, GLOBAL_MATRIZ_DAT_h);


	vTamMuestraDAT = TAMANO_MINIMO_MUESTRA_SEGUNDOS * multiplicadorTiempo;
	vTamMuestraMatriz = vTamMuestraDAT - LEN_FREC_COM;
	huella * listaHuellas = (huella *)malloc((vlenMatriz/vTamMuestraDAT + 1) * sizeof(huella));
	vIndMatriz = 655200;
	clock_t startC = clock();
	while(vIndMatriz < vlenMatriz - vTamMuestraDAT){

		//printf("IND %d - %d , %d , %d\n", vIndMatriz,vlenMatriz - vTamMuestraMatriz,vTamMuestraMatriz, multiplicadorTiempo);
		clock_t startComparador = clock();
		/*if(listaHuellas[vCantidadHuellas].matriz == NULL)
			listaHuellas[vCantidadHuellas].matriz = (entero128 *)malloc(vTamMuestraMatriz*sizeof(entero128));
*/
		listaHuellas[vCantidadHuellas].tamano = vTamMuestraDAT;
		listaHuellas[vCantidadHuellas].inicio = vIndMatriz;

		const int vLimite = TAM_ESCALON*PORCENTAJE_BARRA_MINIMO_PERMITIDO*PUNTOS_ANALISIS_PRIMER_FILTRO;
		int vAvance = vTamMuestraDAT/PUNTOS_ANALISIS_PRIMER_FILTRO;
		if(vAvance == 0)
			vAvance = 1;
		revisarHuella(vIndMatriz, vlenMatriz - 12000, GLOBAL_MATRIZ_DAT_h, &listaHuellas[vCantidadHuellas],vLimite,vAvance,PORCENTAJE_MINIMO_COINCIDENCIAS_PERMITIDO, GLOBAL_ARREGLO_HAMMING_h);

		//printf("AFUERA -- Cantidad Coincidencias  %d \n", listaHuellas[vCantidadHuellas].cantidadCoincidencias);
		//revisarHuella(vIndMatriz , vlenMatriz, vMatrizOUT, &listaHuellas[vCantidadHuellas]);
		if(listaHuellas[vCantidadHuellas].cantidadCoincidencias > 0){
			printf("***************************************************************\nTiempo huella %d/%d =  %f\n",vIndMatriz/multiplicadorTiempo,vlenMatriz/multiplicadorTiempo, ((double)clock() - startComparador)/CLOCKS_PER_SEC);
			printf("Cantidad Coincidencias  %d \n", listaHuellas[vCantidadHuellas].cantidadCoincidencias);
			float tSeg = listaHuellas[vCantidadHuellas].inicio/multiplicadorTiempo;
			float H_hora = tSeg/3600.0f;
			float H_mins = (H_hora - (int)H_hora)*60;
			float H_segs = (H_mins - (int)H_mins)*60;

			for(vIndCoincidencias = 0; vIndCoincidencias < listaHuellas[vCantidadHuellas].cantidadCoincidencias; vIndCoincidencias ++){
				float tSegC = listaHuellas[vCantidadHuellas].coincidencias[vIndCoincidencias].indice/multiplicadorTiempo;
				float C_hora = tSegC/3600.0f;
				float C_mins = (C_hora - (int)C_hora)*60;
				float C_segs = (C_mins - (int)C_mins)*60;
				printf("Huella %d/%d\t%d\t%d -> %d\t%d:%d:%d\t%d -> %d\t%d:%d:%d\t%f\n",vCantidadHuellas,vlenMatriz/vTamMuestraMatriz ,vIndCoincidencias,
				listaHuellas[vCantidadHuellas].inicio,
				listaHuellas[vCantidadHuellas].inicio/multiplicadorTiempo,
				(int)H_hora,(int)H_mins,(int)H_segs,
				listaHuellas[vCantidadHuellas].coincidencias[vIndCoincidencias].indice,
				listaHuellas[vCantidadHuellas].coincidencias[vIndCoincidencias].indice/multiplicadorTiempo,
				(int)C_hora,(int)C_mins,(int)C_segs,
				listaHuellas[vCantidadHuellas].coincidencias[vIndCoincidencias].porcentaje);
			}
			vCantidadHuellas++;
			//HACER CRECER LA MUESTRA
		}else{
			int indMatriz;
			for(indMatriz = 0; indMatriz < vTamMuestraMatriz;  indMatriz ++){
				if(!GLOBAL_MATRIZ_DAT_h[vIndMatriz].vacio)
					GLOBAL_MATRIZ_DAT_h[vIndMatriz + indMatriz].bloqueado = 0;
			}
		}
		//break;


		vIndMatriz += vTamMuestraDAT;

	}



	int i;
	for(i = 0 ; i < vCantidadHuellas;i++ ){
		printf("Indice %d Coincidencias %d\t", i, listaHuellas[i].cantidadCoincidencias);
	}
	printf("\nTAM %d %d\n", lenBuffer, vCantidadHuellas);
}

/*
 *fsize: indica el tamaño en bytes de un archivo.
 * IN: nombre del archivo (const char *)
 * OUT: tamaño en bytes del archivo (int)
 */
int fsize(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1;
}

/*
 leerSuperDat: lee las primeras cuatro frecuencias de un archivo DAT
 IN:
 ** pPathDatIn-> ruta del archivo DAT (char *)
 ** lenBuffer-> tamaño de cada frecuencia (int)
 ** sliceDAT-> variable de salida, matriz donde se copian la informacion de cada frecuencia (unsigned short **)
*/
void leerSuperDat(const char* pPathDatIn, int lenBuffer,unsigned short **  sliceDAT){
    FILE *vArchivo;
    int vIndFrecuency;
    vArchivo = fopen(pPathDatIn,"rb");
    for(vIndFrecuency = 0;vIndFrecuency<NUM_FREC_COM;vIndFrecuency++){
        fread(sliceDAT[vIndFrecuency],NUM_BYTES_SHORT,lenBuffer,vArchivo);
    }
    fclose(vArchivo);
}

/*
 discretizarDAT: lee las primeras cuatro frecuencias de un archivo DAT
 IN:
 ** pLenBuffer-> tamaño de cada frecuencia (int)
 ** pSliceDAT->  informacion de cada frecuencia del dat(unsigned short **)
 ** pMatrizOUT-> variable de salida, matriz donde se almacenan los escalones (unsigned int **)
*/
void discretizarDAT(const int pLenBuffer, unsigned short **pSliceDAT, entero128 *pMatrizOUT){
    const int vLenMatriz = pLenBuffer - LEN_FREC_COM;
    unsigned int vIndCantFrecu = 0;
    int vElemIniEscalon = 0;
    while(vElemIniEscalon < vLenMatriz){
		//printf("ANTES DEL FOR %d - %d\n", vElemIniEscalon , vLenMatriz);
        //unsigned int vEscalon[NUM_FREC_COM] = {0};
        for(vIndCantFrecu = 0; vIndCantFrecu < NUM_FREC_COM;vIndCantFrecu ++){
            unsigned short vArray32Frecu[LEN_FREC_COM] = {0};
            unsigned short vAverageBlock = 0;
            unsigned int vSumValuBlock = 0;
            int vIndFrecuencyBlock = 0;
            short vInd32Block = 0;
            for(vInd32Block = 0; vInd32Block < LEN_FREC_COM;vInd32Block ++){
                vArray32Frecu[vInd32Block]=pSliceDAT[vIndCantFrecu][vElemIniEscalon + vInd32Block];
                vSumValuBlock += vArray32Frecu[vInd32Block];
            }
            //Discretizar los 32 valores
            float value = (((float)vSumValuBlock)/((float)LEN_FREC_COM))*CONST_MEAN[vIndCantFrecu];
            vAverageBlock = (short)value;

            pMatrizOUT[vElemIniEscalon].nums[vIndCantFrecu] = 0;
            for(vIndFrecuencyBlock = 0; vIndFrecuencyBlock < LEN_FREC_COM;vIndFrecuencyBlock++){
                if(vArray32Frecu[vIndFrecuencyBlock] > vAverageBlock){
                    pMatrizOUT[vElemIniEscalon].nums[vIndCantFrecu] <<= 1;
                    pMatrizOUT[vElemIniEscalon].nums[vIndCantFrecu] |= 0b1;
                }
                else
                    pMatrizOUT[vElemIniEscalon].nums[vIndCantFrecu] <<= 1;
            }
        }

        //printf("ANTES DEL MEM %d - %d\n", vElemIniEscalon , vLenMatriz);

        //memcpy(pMatrizOUT[vElemIniEscalon].nums,vEscalon,sizeof(vEscalon));
        //printf("%d\t",pSliceDAT[0][vElemIniEscalon] );
        if(pSliceDAT[0][vElemIniEscalon] == 65535 || pSliceDAT[0][vElemIniEscalon] == 15000|| pSliceDAT[0][vElemIniEscalon] == 39064)
			pMatrizOUT[vElemIniEscalon].vacio = 1;
		else
			pMatrizOUT[vElemIniEscalon].vacio = 0;
		pMatrizOUT[vElemIniEscalon].bloqueado = 0;
        vElemIniEscalon ++;
    }
}

/*
 extraerMuestra: Extrae la matriz de escalones de una seccion del DAT
 IN:
 ** pIndiceInicio-> indice de desplazamiento en el dat (int)
 ** pTamano-> tamaño de la muestra (const int)
 ** pSliceDAT-> Dat completo (unsigned short **)
 ** pMatrizOUT-> variable de salida, matriz donde se almacenan los escalones (entero128 *)

void extraerMuestra(int pIndiceInicio, const int pTamano, entero128 *pMatrizOUT_DAT,  huella *pHuella_Muestra){
	/*unsigned short **sliceMuestra;
    sliceMuestra = (unsigned short **)malloc(NUM_FREC_COM*sizeof(short*));
    int vIndFrecuency;
    for(vIndFrecuency = 0;vIndFrecuency<NUM_FREC_COM;vIndFrecuency++){
        sliceMuestra[vIndFrecuency] = (unsigned short *)malloc(pTamano*sizeof(short));
        memcpy(sliceMuestra[vIndFrecuency], pSliceDAT[vIndFrecuency] + pIndiceInicio, pTamano);
    }

    discretizarDAT(pTamano, sliceMuestra,pMatrizOUT);
    for(vIndFrecuency = 0;vIndFrecuency<NUM_FREC_COM;vIndFrecuency++){
        free(sliceMuestra[vIndFrecuency]);
    }
    free(sliceMuestra);

    //printf("Copiamdo  TAM %d,  \n");
	memcpy(pHuella_Muestra->matriz, pMatrizOUT_DAT + pIndiceInicio, pTamano*sizeof(entero128));
	int indMatriz,vIndEscalon, vDiffHamming;
	for(indMatriz = - multiplicadorTiempo; indMatriz < pTamano;  indMatriz ++){
		if(pIndiceInicio + indMatriz >= 0)
			pMatrizOUT_DAT[pIndiceInicio + indMatriz].bloqueado = 1;

	}
}
*/
/*
 revisarHuella: toma una huella y la compara solo en PUNTOS_ANALISIS_PRIMER_FILTRO cantidad de puntos
 IN:
 ** pIndiceInicio-> inicio del analisis
 ** pLenMatriz-> tamano de la matiz de comparación
 ** pMatrizOUT_DAT-> matriz con escalones del DAT (entero128 *)
 ** pHuella_Muestra-> huella en analisis (huella)
*/
/*
 *
	//printf("indices %d %d\n",pIndiceInicio, pLenMatriz);
	//const int vLimite = TAM_ESCALON*PORCENTAJE_BARRA_MINIMO_PERMITIDO;
			//printf("DESBLOQUEADA %d\n",vIndMatriz);
			for(vIndMatriz_Huella = 0; vIndMatriz_Huella < pHuella_Muestra->tamano; vIndMatriz_Huella += vAvance){
				if(!pMatrizOUT_DAT[vIndMatriz + vIndMatriz_Huella].bloqueado && !pMatrizOUT_DAT[vIndMatriz + vIndMatriz_Huella].vacio && vIndMatriz + vIndMatriz_Huella !=pHuella_Muestra->inicio + vIndMatriz_Huella ){

					vDiffHamming = 0;
					for(vIndEscalon = 0; vIndEscalon < NUM_INT_ESCALON; vIndEscalon ++){
						vDiffHamming += __builtin_popcount(pMatrizOUT_DAT[vIndMatriz + vIndMatriz_Huella].nums[vIndEscalon] ^ pMatrizOUT_DAT[pHuella_Muestra->inicio + vIndMatriz_Huella].nums[vIndEscalon]);
					}


					if(vDiffHamming < vLimite){
						vCantidadValidas ++;
						printf("%d %d\n",vIndMatriz + vIndMatriz_Huella , pHuella_Muestra->inicio + vIndMatriz_Huella );
					}

				}else{
					//printf("BLOQUEADA %d\n",vIndMatriz);
				}
			}
			vPorcentajeSimilitud = (float)vCantidadValidas/PUNTOS_ANALISIS_PRIMER_FILTRO;
			if(vPorcentajeSimilitud >= PORCENTAJE_MINIMO_COINCIDENCIAS_PERMITIDO){
				printf("Limite %f  Porcentaje  %f   validas %d  %d\n", PORCENTAJE_MINIMO_COINCIDENCIAS_PERMITIDO,vPorcentajeSimilitud ,vCantidadValidas, pHuella_Muestra->cantidadCoincidencias);
				pHuella_Muestra->coincidencias[pHuella_Muestra->cantidadCoincidencias].indice = vIndMatriz;
				pHuella_Muestra->coincidencias[pHuella_Muestra->cantidadCoincidencias].porcentaje = vPorcentajeSimilitud;
				pHuella_Muestra->cantidadCoincidencias = (pHuella_Muestra->cantidadCoincidencias + 1)%CANTIDAD_MAXIMA_RECONOCIMIENTOS;
				vIndMatriz += TAMANO_MINIMO_MUESTRA_SEGUNDOS*multiplicadorTiempo;
			}
			*/


void obtenerHamming(int pIndiceInicio, const int pLenMatriz, entero128 *pMatrizOUT_DAT, int pTam_Muestra, int pIni_Muestra, unsigned int * pArregloHamming,const int pAvance){
	int vIndMatriz, vIndMatriz_Huella, vIndEscalon;
	int vDiffHamming;
	for(vIndMatriz =pIndiceInicio; vIndMatriz < pLenMatriz - pTam_Muestra; vIndMatriz ++){
		vDiffHamming = 0;
		for(vIndMatriz_Huella = 0; vIndMatriz_Huella < pTam_Muestra; vIndMatriz_Huella += pAvance){

			for(vIndEscalon = 0; vIndEscalon < NUM_INT_ESCALON; vIndEscalon ++){
				vDiffHamming += __builtin_popcount(pMatrizOUT_DAT[vIndMatriz + vIndMatriz_Huella].nums[vIndEscalon] ^ pMatrizOUT_DAT[pIni_Muestra + vIndMatriz_Huella].nums[vIndEscalon]);
			}

		}
		pArregloHamming[vIndMatriz] = vDiffHamming;

		//break;
	}

}

void revisarHuella(int pIndiceInicio, const int pLenMatriz, entero128 *pMatrizOUT_DAT, huella * pHuella_Muestra, const int pLimite, const int pAvance, const float pPORCENTAJE,unsigned int * pArregloHamming){
/*
	const int vLimite = TAM_ESCALON*PORCENTAJE_BARRA_MINIMO_PERMITIDO*PUNTOS_ANALISIS_PRIMER_FILTRO;
	int vAvance = pHuella_Muestra->tamano/PUNTOS_ANALISIS_PRIMER_FILTRO;
	if(vAvance == 0)
		vAvance = 1;
*/
	int vIndMatriz, vIndEscalon;
	pHuella_Muestra->cantidadCoincidencias = 0;

	float vPorcentajeSimilitud;


	obtenerHamming(pIndiceInicio,pLenMatriz,pMatrizOUT_DAT, pHuella_Muestra->tamano, pHuella_Muestra->inicio, pArregloHamming, pAvance);
	//printf("Suma Hamming %d \t Limite %d \n",vDiffHamming, vLimite );
	for(vIndMatriz =pIndiceInicio; vIndMatriz < pLenMatriz - pHuella_Muestra->tamano; vIndMatriz ++){
		if((vIndMatriz < pHuella_Muestra->inicio - RESTRICCION_CERCANIA_SEGUNDOS*multiplicadorTiempo ||
		vIndMatriz > pHuella_Muestra->inicio + RESTRICCION_CERCANIA_SEGUNDOS*multiplicadorTiempo)
		&& vIndMatriz  !=pHuella_Muestra->inicio ){
			if(!pMatrizOUT_DAT[vIndMatriz].bloqueado && !pMatrizOUT_DAT[vIndMatriz].vacio ){
				if(pArregloHamming[vIndMatriz] < pLimite){
					vPorcentajeSimilitud = 1 - (float)pArregloHamming[vIndMatriz]/(float)(TAM_ESCALON*PUNTOS_ANALISIS_PRIMER_FILTRO);
					if( vPorcentajeSimilitud > pPORCENTAJE){
						//printf("Limite %f  Porcentaje  %f   validas %d  %d\n", PORCENTAJE_MINIMO_COINCIDENCIAS_PERMITIDO,vPorcentajeSimilitud ,vCantidadValidas, pHuella_Muestra->cantidadCoincidencias);
						pHuella_Muestra->coincidencias[pHuella_Muestra->cantidadCoincidencias].indice = vIndMatriz;
						pHuella_Muestra->coincidencias[pHuella_Muestra->cantidadCoincidencias].porcentaje = vPorcentajeSimilitud;
						pHuella_Muestra->cantidadCoincidencias = (pHuella_Muestra->cantidadCoincidencias + 1)%CANTIDAD_MAXIMA_RECONOCIMIENTOS;
						//vIndMatriz += TAMANO_MINIMO_MUESTRA_SEGUNDOS*multiplicadorTiempo;
					}
				}
			}
		}
	}
	int a = 0;
	//printf("Coincidencias encontradas = %d \n", pHuella_Muestra->cantidadCoincidencias);
}


void incrementarHuella(const int pLenMatriz, entero128 *pMatrizOUT_DAT, huella * pHuella_Muestra,unsigned int * pArregloHamming){
	//hola soy diego

	int a = 4;
}
