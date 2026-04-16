#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include "../../common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>

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

  double G, R_T, R_0, m_A, d, h, m_T, m_L, v_0, rho_O, lambda, Cx, d_TL;         // accélération gravitationnelle, masse, longueur, fréquence angulaire, rayon, coefficient de frottement


std::valarray<double> y ; // tableau de variables d'état du système


  double tf;          // Temps final
  double dt;      // Intervalle de temps
  double numBodies;    // Nombre de corps dans la simulation
  int nsteps_per; // Nombre de pas de temps par période d'excitation

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
      double emec = Emec(y); // TODO: Evaluer l'energie mecanique
      double pnc = Pnonc(y); // TODO: Evaluer la puissance des forces non conservatives
      *outputFile << t << " " << y << " " << emec << " " << pnc << endl;
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

  double mass(size_t i) const
    {
        if (i == 0) return m_A;
        if (i == 1) return m_T;
        return m_L;
    }

  valarray<double> vectxy(size_t i, size_t j, const valarray<double>& y) const
    {
        valarray<double> r(2);
        r[0] = y[ix(i)] - y[ix(j)];
        r[1] = y[iy(i)] - y[iy(j)];
        return r;
    }

  valarray<double> vectv(size_t i, const valarray<double>& y) const
    {
        valarray<double> v(2);
        v[0] = y[ivx(i)];
        v[1] = y[ivy(i)];
        return v;
    }

  valarray<double> vectvxy(size_t i, size_t j, const valarray<double>& y) const
    {
        valarray<double> v(2);
        v[0] = y[ivx(i)] - y[ivx(j)];
        v[1] = y[ivy(i)] - y[ivy(j)];
        return v;
    }

  double normev(const valarray<double>& v) const
    {
        return sqrt(pow(v[0], 2) + pow(v[1], 2));
    }

  double distance(std::size_t i, std::size_t j, const std::valarray<double>& y)
    {
        return sqrt(pow(vectxy(i, j, y)[0], 2) + pow(vectxy(i, j, y)[1], 2));
    }

  double vitesse(std::size_t i, std::size_t j, const std::valarray<double>& y) 
    {
        return sqrt(pow(vectvxy(i, j, y)[0], 2) + pow(vectvxy(i, j, y)[1], 2));
    }


  double Emec(const std::valarray<double>& y)
    {
        return (1/2.0)*m_A*(pow(vectv(1,y)[0],2)+pow(vectv(1,y)[1],2)) + (1/2.0)*m_T*(pow(vectv(2,y)[0],2)+pow(vectv(2,y)[1],2)) + (1/2.0)*m_L*(pow(vectv(3,y)[0],2)+pow(vectv(3,y)[1],2)) - G*m_A*m_T/(distance(1,2,y)) - G*m_A*m_L/(distance(1,3,y)) - G*m_T*m_L/(distance(2,3,y));
    }


  std::valarray<double> qte(const std::valarray<double>& y)
    {
        return m_A*(vectv(1,y)) + m_T*(vectv(2,y)) + m_L*(vectv(3,y)); // Quantité de mouvement du système
    }


  // TODO definir la puissance des forces non conservatives
  double Pnonc( const std::valarray<double>& y  )
  {
      return -(1/2)*rho_O*exp((R_T-distance(1,2,y))/lambda)*S*Cx*pow(vitesse(1,2,y),3) -(1/2)*rho_O*exp((R_T-distance(2,3,y))/lambda)*S*Cx*pow(vitesse(2,3,y),3) ;
  }

  double forcegravitation(std::size_t i, std::size_t j, const std::valarray<double>& y , double a)
  {
      return -G*m_A*a*vectxy(i,j,y)/(pow(distance(i,j,y),3));
  }
  // TODO écrire la fonction pour l'acceleration
  std::valarray<double> acc = [(forcegravitation(1,2,y,m_T) + forcegravitation(1,3,y,m_L))/m_A , (forcegravitation(2,1,y,m_A) + forcegravitation(2,3,y,m_L))/m_T , (forcegravitation(3,1,y,m_A) + forcegravitation(3,2,y,m_T))/m_L];
  

  // TODO implementer le schéma Velocity Verlet pour une accélération dependante du theta, thetadot et t.
  void step()
  {
    double yold = y;
    
    double demiv = thetadot + compute_acc(thetaold, thetadot, t)/2;
     
    theta += dt*thetadot + (dt*dt/2)*compute_acc(theta, thetadot, t);
    
    thetadot+= (dt/2)*(compute_acc(theta, demiv, t) + compute_acc(thetaold, demiv, t));
    
    t += dt;
  }

  void std::valarray<double> rk4Step(double step, const std::valarray<double>& y)
  {

    std::valarray<double> k1 = (y);
    std::valarray<double> k2 = (y + 0.5*k1);
    std::valarray<double> k3 = (y + 0.5*k2);
    std::valarray<double> k4 = (y + k3);
    
    return y + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
  }

void CheckCollision(const std::valarray<double>& y)
{

}


public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe

      m_A     = configFile.get<double>("m_A", m_A);	        
      d     = configFile.get<double>("d", d);         
      R_0    = configFile.get<double>("R_0", R_0);        
      v_0    = configFile.get<double>("v_0", v_0);         
      h = configFile.get<double>("h", h); 
      G     = configFile.get<double>("G", G);        
      R_T = configFile.get<double>("R_T", R_T); 
      m_T    = configFile.get<double>("m_T", m_T);    
      m_L    = configFile.get<double>("m_L", m_L);    
      Cx = configFile.get<double>("Cx", Cx); 
      lambda = configFile.get<double>("lambda", lambda);
      d_TL = configFile.get<double>("d_TL", d_TL);

      x_1 = configFile.get<double>("x_1", x_1);
      x_2 = configFile.get<double>("x_2", x_2);
      x_3 = configFile.get<double>("x_3", x_3);
      y_1 = configFile.get<double>("y_1", y_1);
      y_2 = configFile.get<double>("y_2", y_2);
      y_3 = configFile.get<double>("y_3", y_3);    
      v_x1 = configFile.get<double>("v_x1", v_x1);
      v_x2 = configFile.get<double>("v_x2", v_x2);
      v_x3 = configFile.get<double>("v_x3", v_x3);
      v_y1 = configFile.get<double>("v_y1", v_y1);
      v_y2 = configFile.get<double>("v_y2", v_y2);
      v_y3 = configFile.get<double>("v_y3", v_y3);


      N_excit  = configFile.get<int>("N");            // number of periods of excitation
      nsteps_per= configFile.get<int>("nsteps");        // number of time step per period
      sampling = configFile.get<unsigned int>("sampling",sampling); // lire le nombre de pas de temps entre chaque ecriture des diagnostics

      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output").c_str());
      outputFile->precision(15);
      if(N_excit>0){
        tf = N_excit*(2*pi/Omega);
        dt   = (2*pi/Omega)/nsteps_per;
      }
      else{
        dt = tf/nsteps_per;
      }
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

