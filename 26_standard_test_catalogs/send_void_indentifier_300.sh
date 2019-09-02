#!/bin/bash

#SBATCH --job-name=Vcirmax300           #Nombre del job
#SBATCH -p medium			#Cola a usar, Default=short (Ver colas y límites en /hpcfs/shared/README/partitions.txt)
#SBATCH -N 1				#Nodos requeridos, Default=1
#SBATCH -n 1				#Tasks paralelos, recomendado para MPI, Default=1
#SBATCH --cpus-per-task=1		#Cores requeridos por task, recomendado para multi-thread, Default=1
#SBATCH --mem-per-cpu=2048		#Memoria en Mb por CPU, Default=2048
#SBATCH --time=159:00:00		#Tiempo máximo de corrida, Default=2 horas. Formato HH:MM:SS
#SBATCH --mail-user=fl.gomez10@uniandes.edu.co
#SBATCH --mail-type=ALL			
#SBATCH -o Vcircmax_300_.o%j.LOG	#Nombre de archivo de salida

host=`/bin/hostname`
date=`/bin/date`
echo "BetaSkeletonVoidFinder en HPC Magnus"
echo "Corri en la maquina: "$host
echo "Corri el: "$date

module load anaconda/python3.7

make 
