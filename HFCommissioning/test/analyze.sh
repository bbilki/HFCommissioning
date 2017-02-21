#!/bin/sh

if [ ! -d Data ]
then
	mkdir Data
fi

if [ ! -d Histos ]
then
	mkdir Histos
fi

if [ ! -d Plots ]
then
	mkdir Plots
fi

if [ ! -d NTuples ]
then
	mkdir NTuples
fi

cmsRun hfcommissioning_cfg.py $1 $2
wait
mv N_$1.root NTuples
wait

cd Files
./runPlotter.sh $1 $2
wait
cd ..

cd Plots
tar -cf plots.tar $1
wait
gzip plots.tar
wait
scp -P 53222 plots.tar.gz hfSX5@feynman.physics.uiowa.edu:/var/www/html/HFCommissioning2017
wait
ssh -p 53222 hfSX5@feynman.physics.uiowa.edu "cd /var/www/html/HFCommissioning2017 ; tar -zxf /var/www/html/HFCommissioning2017/plots.tar.gz ; wait ; rm plots.tar.gz ; ./makePlotList.sh"
wait
rm plots.tar.gz
cd -

