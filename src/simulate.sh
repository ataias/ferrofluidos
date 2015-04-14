#!/usr/bin/env bash

#Faz simulação e gera vídeo para centro do campo magnético em (a,b)
#0 <= a <= 1, mesmo para b
simul () {
a=$1
b=$2
FOLDERTOSAVE=$3
VIDEONAME=$4 #no extension
rm -rf $FOLDERTOSAVE
mkdir -p $FOLDERTOSAVE

#Parâmetros
n=22; #tamanho da malha
t=1.0;
Re=10;
divFactor=1.25;
chi=0.5;
Cpm=0.8;
save=1; #para salvar o arquivo de simulação, deve ser setado em 1

FILENAME=N$((n - 2)).dat

step=1; #para plotar png

#simulação
julia frames.jl $n $t $Re $divFactor $chi $Cpm $save $a $b > out$((n-2)).txt
rm -rf png
mkdir -p png && cd png
mv -v ../$FILENAME .
../vectorField.py $((n-2)) $t $step $chi $Cpm $Re
mv -v $FILENAME $FOLDERTOSAVE/
cd ..
mv -v out$((n-2)).txt $FOLDERTOSAVE/
mv -v png $FOLDERTOSAVE/png
ffmpeg -i $FOLDERTOSAVE/png/%4d.png -c:v libx264 -vf fps=30 -pix_fmt yuv420p $FOLDERTOSAVE/$VIDEONAME.mp4
}

#simul 0.0 0.5 /Users/ataias/Documents/ferrisimulacao/sim01 sim01
simul 0.5 0.5 /Users/ataias/Documents/ferrisimulacao/sim02 sim02
simul 0.5 1.0 /Users/ataias/Documents/ferrisimulacao/sim03 sim03
simul 0.0 0.0 /Users/ataias/Documents/ferrisimulacao/sim04 sim04