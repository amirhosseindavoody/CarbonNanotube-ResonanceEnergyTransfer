#!/bin/bash
clear
cwd=$(pwd)

#n_ch=07
#m_ch=05
#nkg=1001
#nr=0200
#E_th=1.5
#Kcm_max=1.5
#flg_dielectric=1
#i_sub=1
#kappa=6.055

#echo "Chirality is ($n_ch , $m_ch)"
#echo "Dielectric flag is" $flg_dielectric
#echo nkg = $nkg
#echo E_th= $E_th
#echo nr= $nr
#echo Kcm_max= $Kcm_max
#echo i_sub= $i_sub
#echo kappa= $kappa

cd "/cygdrive/c/Users/Amirhossein/Google Drive/Research/Exciton/Data/Environmental Effect"
#targefolder="CNT($n_ch,$m_ch)-nkg($nkg)-nr($nr)-E_th($E_th)-Kcm_max($Kcm_max)-i_sub($i_sub)-kappa($kappa)"

#answer='y'
#if [ -d "$targefolder" ] ; then
#	echo "Directory already exists!"
#	echo "Do you want to redo the simulation? (y/n)"
#	read answer
#else
#	mkdir $targefolder
#fi

#if [ "$answer" = "y" ] ; then
echo "Simulation started in background!"
#	rm -rf $targefolder
#	mkdir $targefolder
#	cd "$targefolder"
$cwd/main.exe
#else
#	echo "Simulation NOT started!"
#fi
