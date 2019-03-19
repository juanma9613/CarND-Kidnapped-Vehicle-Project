/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"
#include <limits>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
// why static???
static std::default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 10;  // TODO: Set the number of particles
 
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  // recall : std::vector<Particle> particles; it's a vector of particles
  for (int i=0; i<num_particles; i++)  
  {
  Particle temporal_particle;
  temporal_particle.id = i;
  temporal_particle.x = dist_x(gen);
  temporal_particle.y = dist_y(gen);
  temporal_particle.theta = dist_theta(gen);
  temporal_particle.weight=1;
  // adding gaussian noise  
  weights.push_back(1);  
  particles.push_back(temporal_particle);

  }
  
  
  //std::cout<<particles[0].weight<<"particle 0 weight w\n";
  is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  double new_x;
  double new_y;
  double new_theta;
  
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);
  for (int i=0;i<num_particles;i++)
  {
    if ( fabs(yaw_rate) < 0.01)
    {
      new_x = particles[i].x + (velocity * delta_t * cos(particles[i].theta));
      new_y = particles[i].y + (velocity * delta_t * sin(particles[i].theta));
      new_theta = particles[i].theta;
     }
    else
    {
      new_theta = particles[i].theta + delta_t*yaw_rate;
      new_x = particles[i].x + velocity/yaw_rate*(sin(new_theta) -sin(particles[i].theta));
      new_y = particles[i].y + velocity/yaw_rate*(-cos(new_theta) +cos(particles[i].theta));     
    }
  particles[i].x = dist_x(gen) +new_x;
  particles[i].y = dist_y(gen) +new_y;
  particles[i].theta = dist_theta(gen)+new_theta;    
  //std::cout<<"Predicting new value";
  //std::cout<<particles[i].x<<"PREDICTED particles x"<<i<<"\n";
  //std::cout<<particles[i].y<<"PREDICTED particles y"<<i<<"\n";
  //std::cout<<particles[i].theta<<"PREDICTED particles theta"<<i<<"\n";
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) 
                                  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase..
   */
{     // for each observation called o passed by reference
  
  	 // TODO in the future: let's think what would happen if the predicted vector is null?
      for(LandmarkObs& o: observations)
      {
      double min = std::numeric_limits<double>::max();    
            for(LandmarkObs& p: predicted)
            {
              double d = dist(o.x, o.y, p.x, p.y);
                    if(d < min)
                    {	
                  min = d;
                      // changing the id of the observation to match the id of the landmark predicted
                  o.id = p.id;
                    }
            }    
      } 


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) 
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    for (int i=0; i<num_particles; i++)  
  {
      // get all observations in the range of the particle
          //check range, get close landmarks
    vector<LandmarkObs> predictions;    
      for(unsigned int landmark = 0 ;landmark<map_landmarks.landmark_list.size();landmark++) 
      {
              double d = dist(map_landmarks.landmark_list[landmark].x_f, map_landmarks.landmark_list[landmark].y_f, particles[i].x, particles[i].y); 
              if(d < sensor_range)
              {
                LandmarkObs new_obs = {};
                new_obs.id = map_landmarks.landmark_list[landmark].id_i;
                new_obs.x  = map_landmarks.landmark_list[landmark].x_f;
                new_obs.y  = map_landmarks.landmark_list[landmark].y_f;
                predictions.push_back(new_obs);
         		}
      }
      
      vector<LandmarkObs> observations_trans;
      // transform each observation to the map coordinate system
          for (unsigned int j=0;j<observations.size();j++)
          {
            LandmarkObs obs_trans;
              obs_trans.id = observations[j].id;
              obs_trans.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) -    (sin(particles[i].theta) * observations[j].y);
              obs_trans.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
            observations_trans.push_back(obs_trans);
          }
  dataAssociation(predictions,observations_trans);
      // after this step, each observation should have x,y in world position and should also have a pair in the observations
      // Now the important thing is to update the weights
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;
    //compute weights
    double weight=1;
    double std_2_pi = 2.0*M_PI*std_landmark[0]*std_landmark[1];
    double std_x_2  = 2.0*std_landmark[0]*std_landmark[0];
    double std_y_2 = 2.0*std_landmark[1]*std_landmark[1];
    //std::cout<<"values of sd"<< std_2_pi<<"\t"<<std_x_2<<"\t"<<std_y_2;
     // iteration the vector obs_trans where each is an struct of the landmark type 
    for(LandmarkObs& obs : observations_trans)
    {
      // remeber id = 1 is the landmarklist 0 for the map
      Map::single_landmark_s mark = map_landmarks.landmark_list[obs.id -1];
      double e1 = std::pow(obs.x - mark.x_f, 2);
      double e2 = std::pow(obs.y - mark.y_f, 2);
     // std::cout<<"\n values of e1 e2 e"<< e1<<"\t"<<e2<<"\t";
      double e = (e1/std_x_2 + e2/std_y_2);

      double ee = exp(-e);
      double w = ee/std_2_pi;
    //  std::cout<<"\n values of e ee w"<< e<<"\t"<<ee<<"\t"<<w;
      //prod of all weights
      weight *= w;
      //record association
      associations.push_back(obs.id);
      sense_x.push_back(obs.x);
      sense_y.push_back(obs.y); 
    }
    particles[i].weight = weight;

      
    //insert into weight vector
    weights[i]= weight;
    //update particle's associations
SetAssociations(particles[i], associations, sense_x, sense_y); 
  }
    // for (int gh = 0;gh<num_particles;gh++)
    //{std::cout<<particles[gh].weight<<"particles weights";
    //}
  //std::cout<<"UPDATING W\n";
}

void ParticleFilter::resample() {
  
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
     // std::random_device rd;
    // std::mt19937 gen(rd()); not necessary if the static is used
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
   vector<Particle> n_particles;
    //resample  
  for(int i=0; i < num_particles; i++){      
    const int index = distribution(gen);  
    n_particles.push_back(particles[index]);
  }
  //clear Particle vectors
  //particles.clear();
  
  //new set of particles
  particles = n_particles;
  //std::cout<<"resampling particles w\n";
   

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}