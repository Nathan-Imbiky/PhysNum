#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include "../../common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>
#include <valarray>


using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private:
    // Existing private members of Engine...
  const double pi=3.1415926535897932384626433832795028841971e0;

  // definition des variables

  double G, mA, d, r0, v0, h, mT, mL, rho_0, R_T, lambda, Cx, dTL;         // accélération gravitationnelle, masse, longueur, fréquence angulaire, rayon, coefficient de frottement

  bool adaptative = true;
  std::valarray<double> y;

  double t;  // Temps courant pas de temps
  double tf;          // Temps final
  double dt;      // Intervalle de temps
  int nsteps_per; // Nombre de pas de temps par période d'excitation
  int numBodies;

  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      double emec = Emec(y, t); // TODO: Evaluer l'energie mecanique
      double pnc = Pnonc(y, t); // TODO: Evaluer la puissance des forces non conservatives
      *outputFile << t << " " << y[ix(0)] << " " << y[iy(0)] << " " << y[ivx(0)] << "" << y[ivy(0)] <<"" << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }
  
  std::size_t ix(std::size_t i) const { return 2 * i; } 
  std::size_t iy(std::size_t i) const { return 2 * i + 1; } 
  std::size_t ivx(std::size_t i) const { return 2 * numBodies + 2 * i; } 
  std::size_t ivy(std::size_t i) const { return 2 * numBodies + 2 * i + 1; }
  
  // TODO definir l'énergie mechanique
  double Emec(std::valarray<double> const& y, double t_)
  {
      double K =(1/2.0)*(mA*(y[ivx(0)]*y[ix(0)]) + mT*(y[ivx(1)]*y[ivx(1)]) + mL*(y[ivx(2)]*y[ivx(2)]));
      double U = -G*( mA*mT/( norm(0, 1) ) + mA*mL/( norm(0, 2) ) + mL*mT/( norm(1, 2) ));
      return  K + U;
  }
  
	  double norm(size_t i, size_t j) ////distance entre les astres i et j
	  {
		  if(i<y.size()/6 && j<y.size()/6)
		  {
			return sqrt((y[ix(i)]-y[ix(j)])*(y[ix(i)]-y[ix(j)]) + (y[iy(i)]-y[iy(j)])*(y[iy(i)]-y[iy(j)]));
		}
		else
		{
			return 0.0;
		}
	  }
  double Pnonc(std::valarray<double> const& y, double t_)
  {
	  return ;
  }
  
  std::valarray<double> momentum(std::valarray<double> const& y, double t_)
  {
	  return mA*y[std::slice(ivx(0), 2, 1)] + mT*y[std::slice(ivx(1), 2, 1)] + mL*y[std::slice(ivx(2), 2, 1)];
  }

checkCollisions, printOut, rk4step, step, run, Pnonc, jspquoi d'autre mais au moins tt ça qui reste 

  // TODO definir la puissance des forces non conservatives
  double Pnonc(std::valarray<double> const& y, double t_)
  {
      
      return 0.;
  }
  
  double mass(size_t i) ////helper associant la masse d'un astre à son numéro
	  {
		  
        switch (i)
		{
		case 0:
			return mA;
		case 1:
			return mT;
		case 2:
			return mL;
		default:
			return 0.0;
		}
	  }
  
  double gravx(size_t i, size_t j) ////helper pour l'expression de la force de gravitation sur le corps i par le corps j selon l'axe x
	  {
		  if(i<y.size()/6 && j<y.size()/6 && i!=j)
		  {
			return -G*mass(i)*mass(j)*(y[ix(i)]-y[ix(j)])/(norm(i, j)*norm(i, j)*norm(i, j));
		}
		else
		{
			return 0.0;
		}
	  }
	  
  double gravy(size_t i, size_t j) ////helper pour l'expression de la force de gravitation sur le corps i par le corps j selon l'axe y
	  {
		  if(i<y.size()/6 && j<y.size()/6 && i!=j)
		  {
			return -G*mass(i)*mass(j)*(y[iy(i)]-y[iy(j)])/(norm(i, j)*norm(i, j)*norm(i, j));
		}
		else
		{
			return 0.0;
		}
	  }	  

  // TODO écrire la fonction pour l'acceleration (theta_doubledot)
  std::valarray<double> acc(std::valarray<double> const& y)
  {
      std::valarray<double> f = {y[ivx(0)], y[ivy(0)], y[ivx(1)], y[ivy(1)], y[ivx(2)], y[ivy(2)], gravx(0, 1) + gravx(0, 2), gravy(0, 1) + gravy(0, 2), gravx(1, 0) + gravx(1, 2), gravy(1, 0) + gravy(1, 2), gravx(2, 1) + gravx(2, 0), gravy(2, 1) + gravy(2, 0),};

      return f;
  }
  // TODO implementer le schéma Velocity Verlet pour une accélération dependante du theta, thetadot et t.
  void step()
  {
	
    t += dt;
  }


public:
    // Modified constructor
    Engine(ConfigFile configFile)///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      tf     = configFile.get<double>("tf",tf);	        // t final (overwritten if N_excit >0)
      G     = configFile.get<double>("G", G);         // lire l'acceleration de gravite
      mA     = configFile.get<double>("mA", mA);         // lire la masse
      d = configFile.get<double>("d", d);
      r0 = configFile.get<double>("r0", r0);
      v0 = configFile.get<double>("v0", v0);
      h = configFile.get<double>("h", h);
      mT = configFile.get<double>("mT", mT);
      mL = configFile.get<double>("mL", mL);
      rho_0 = configFile.get<double>("rho_0", rho_0);
      R_T = configFile.get<double>("R_T", R_T);
      lambda = configFile.get<double>("lambda", lambda);
      Cx = configFile.get<double>("Cx", Cx);
      dTL = configFile.get<double>("dTL", dTL);  
      
     double x1= configFile.get<double>("x1", x1);
     double x2= configFile.get<double>("x2", x2);
     double x3= configFile.get<double>("x3", x3);
     double y1= configFile.get<double>("y1", y1);
     double y2= configFile.get<double>("y2", y2);
     double y3= configFile.get<double>("y3", y3);
     double vx1= configFile.get<double>("vx1", vx1);
     double vx2= configFile.get<double>("vx2", vx2);
     double vx3= configFile.get<double>("vx3", vx3);
     double vy1= configFile.get<double>("vy1", vy1);
     double vy2= configFile.get<double>("vy2", vy2);
     double vy3= configFile.get<double>("vy3", vy3);
      
      y = {x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3};

      nsteps_per= configFile.get<int>("nsteps");        // number of time step per period
      sampling = configFile.get<unsigned int>("sampling",sampling); // lire le nombre de pas de temps entre chaque ecriture des diagnostics

      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output").c_str());
      outputFile->precision(15);
      
      dt = tf/nsteps_per;
      numBodies = y.size();
    };


    // Destructeur virtuel
    virtual ~Engine()
    {
      outputFile->close();
      delete outputFile;
    };
      // Simulation complete
    void run()
    {
      t = 0.;
      last = 0;
      printOut(true);

      while( t < tf-0.5*dt )
      {
        step();
        printOut(false);
      }
      printOut(true);
    };
};

// programme
int main(int argc, char* argv[])
{
  // Existing main function implementation
  // ...
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

  Engine* engine;

  // Create an instance of Engine instead of EngineEuler
  engine = new Engine(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}


