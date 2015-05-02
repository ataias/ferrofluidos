#!/usr/bin/env bash

#Faz simulação e gera vídeo para centro do campo magnético em (a,b)
#0 <= a <= 1, mesmo para b
simul () {
Re=$1
a=$2
b=$3
FOLDERTOSAVE=$4
VIDEONAME=$5 #no extension
gamma=$6
rm -rf $FOLDERTOSAVE
mkdir -p $FOLDERTOSAVE

#Parâmetros
n=32; #tamanho da malha
t=1.5;
divFactor=1.25;
chi=0.5;
Cpm=0.8;

save=1; #para salvar o arquivo de simulação, deve ser setado em 1

FILENAME=N$((n - 2)).dat

step=1; #para plotar png

#simulação
julia frames.jl $n $t $Re $divFactor $chi $Cpm $save $a $b $gamma > out$((n-2)).txt
rm -rf png
mkdir -p png && cd png
mv -v ../$FILENAME .
../vectorField.py $((n-2)) $t $step $chi $Cpm $Re $gamma
mv -v $FILENAME $FOLDERTOSAVE/
cd ..
mv -v out$((n-2)).txt $FOLDERTOSAVE/
mv -v png $FOLDERTOSAVE/png
#ffmpeg -i $FOLDERTOSAVE/png/%4d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p $FOLDERTOSAVE/$VIDEONAME.mp4
}

#Simulações sem magnetismo, gamma = 0
gamma=0
#Re=1; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=10; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
Re=20; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=30; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=40; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=50; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=60; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=70; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=80; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=90; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma
#Re=100; simul $Re 0.0 -0.05 ~/Documents/simulacao/Re$Re Re$Re $gamma

gamma=5
#Re=1; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
Re=10; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=20; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=30; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=40; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=50; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=60; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=70; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=80; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=90; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma
#Re=100; simul $Re 0.0 -0.05 ~/Documents/simulacao/mRe$Re Re$Re $gamma