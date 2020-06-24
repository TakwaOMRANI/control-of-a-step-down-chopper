%Période d'échantillonnage
Te=20e-6;
%Période MLI
TMLI=200e-6;
%Pramètres du système
R=0.3;
L=12.5e-3;
Cdc=200e-6;
Rch=50;
Ve=200;
%Caractéristique de la réponse du système en BO à un échelon appliqué à
%Vconv
GBO=Rch/(R+Rch);
wnBO=sqrt((R+Rch)/(L*Cdc*Rch));
KsiBO=(1/(2*wnBO))*(R/L+1/(Cdc*Rch));
p1BO=-KsiBO*wnBO+j*wnBO*sqrt(1-KsiBO^2);
p2BO=-KsiBO*wnBO-j*wnBO*sqrt(1-KsiBO^2);
DBO=100*exp(-pi*KsiBO/sqrt(1-KsiBO^2));
tpBO=pi/(wnBO*sqrt(1-KsiBO^2));
%Matrices du modèle d'état du système
A=[-R/L -1/L;1/Cdc -1/(Cdc*Rch)];
B=[1/L;0];
C=[0 1];
X_0=[0;0];
%Synthèse de la commande par Retour d'Etat
DBF=10;
tpBF=2e-3;
KsiBF=sqrt((log(DBF/100)^2)/(pi^2+log(DBF/100)^2));
wnBF=pi/(tpBF*sqrt(1-KsiBF^2));
p1BF=-KsiBF*wnBF+j*wnBF*sqrt(1-KsiBF^2);
p2BF=-KsiBF*wnBF-j*wnBF*sqrt(1-KsiBF^2);
K=place(A,B,[p1BF;p2BF]);
N=inv(-C*inv(A-B*K)*B);
%Synthèse de l'observateur de Luenberger
O=[C;C*A];
p1OBS=2*p1BO;
p2OBS=2*p2BO;

% p1OBS=-0.7*5*wnBF+j*5*wnBF*sqrt(1-0.7^2);
% p2OBS=-0.7*5*wnBF-j*5*wnBF*sqrt(1-0.7^2);

Lamda1_OBS=p1OBS;
Lamda2_OBS=p2OBS;
LOBS_T=place(A',C',[Lamda1_OBS;Lamda2_OBS]);
LOBS=LOBS_T';
%Observateur de Luenberger discret
A=[-R/L -1/L;1/Cdc -1/(Cdc*Rch)];
B=[1/L;0];
C=[0 1];
Ad=[1-Te*R/L -Te/L;Te/Cdc 1-Te/(Cdc*Rch)];
Bd=[Te/L;0];
Cd=[0 1];
LOBSd=Te*LOBS;