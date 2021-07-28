/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
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
using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  // define normal distributions for sensor noise
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    Particle current_p;
    current_p.id = i;
    current_p.x = dist_x(gen);
    current_p.y = dist_y(gen);
    current_p.theta = dist_theta(gen);
    current_p.weight = 1.0;

    particles.push_back(current_p);
    weights.push_back(current_p.weight);
  }
  is_initialized = true;
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

   for (unsigned int i = 0; i < num_particles; i++) {
     double particle_x = particles[i].x;
     double particle_y = particles[i].y;
     double particle_theta = particles[i].theta;

     double pred_x;
     double pred_y;
     double pred_theta;
      // Instead of a hard check of 0, adding a check for very low value of yaw_rate
     if (fabs(yaw_rate) < 0.00001) {
          pred_x = particle_x + velocity * cos(particle_theta) * delta_t;
          pred_y = particle_y + velocity * sin(particle_theta) * delta_t;
          pred_theta = particle_theta;
     } else {
          pred_x = particle_x + (velocity / yaw_rate) * (sin(particle_theta + (yaw_rate * delta_t)) - sin(particle_theta));
          pred_y = particle_y + (velocity / yaw_rate) * (cos(particle_theta) - cos(particle_theta + (yaw_rate * delta_t)));
          pred_theta = particle_theta + (yaw_rate * delta_t);
     }
     // define normal distributions for sensor noise
      normal_distribution<double> dist_x(pred_x, std_pos[0]);
      normal_distribution<double> dist_y(pred_y, std_pos[1]);
      normal_distribution<double> dist_theta(pred_theta, std_pos[2]);
      // add noise
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
   }  
}
void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

   for (unsigned int i = 0; i < observations.size(); i++) {
      // get current observation object 
      LandmarkObs obser = observations[i];
      // init minimum distance to maximum double value
      double nearest_dist = numeric_limits<double>::max();
      int nearest_landmark_id = -1;
      
      for (unsigned int j = 0; j < predicted.size(); j++) {
        // get current prediction object
        LandmarkObs pre = predicted[j];
        double current_dist = dist(obser.x, obser.y, pre.x, pre.y);
        // find the nearest distance to the current observed landmark
        if (current_dist < nearest_dist) {
              nearest_dist = current_dist;
              nearest_landmark_id = pre.id;
        }
      }
      observations[i].id = nearest_landmark_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
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

   /*This variable is used for normalizing weights of all particles and bring
      them in the range of [0, 1]*/
   double weight_normalizer = 0.0;

   for (unsigned int i = 0; i < num_particles; i++) {
     double particle_x = particles[i].x;
     double particle_y = particles[i].y;
     double particle_theta = particles[i].theta;

     /*Step 1: Transform observations from vehicle coordinates to map coordinates.*/
     vector<LandmarkObs> transformed_observations;
      // create and populate a copy of the list of observations transformed from vehicle coordinates to map coordinates
     for (unsigned int j = 0; j < observations.size(); j++) {
        double transform_x = particle_x + (cos(particle_theta) * observations[j].x) - (sin(particle_theta) * observations[j].y);
        double transform_y = particle_y + (sin(particle_theta) * observations[j].x) + (cos(particle_theta) * observations[j].y);
        transformed_observations.push_back(LandmarkObs {observations[j].id, transform_x, transform_y});
     }

     /*Step 2: Filter map landmarks to only keep those which are in the sensor_range of current particle and push them to predictions vector.*/
     vector<LandmarkObs> predicted_landmarks;
     for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
          Map::single_landmark_s current_landmark = map_landmarks.landmark_list[j];
          if ((fabs((particle_x - current_landmark.x_f)) <= sensor_range) && (fabs((particle_y - current_landmark.y_f)) <= sensor_range)) {
                predicted_landmarks.push_back(LandmarkObs{current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
          }
     }

     /*Step 3: Using nearest neighbor algorithm to get associate observations to predicted landmarks .*/

     dataAssociation(predicted_landmarks, transformed_observations);

    /*Step 4: Calculate the weight of each particle by using Multivariate Gaussian Distribution.*/
    // Reset the weight of particle to 1.0
    particles[i].weight = 1.0;

     /*Calculate the weight of particle by using multivariate Gaussian probability function*/
    for (unsigned int k = 0; k < transformed_observations.size(); k++) {
       double t_obs_x = transformed_observations[k].x;
       double t_obs_y = transformed_observations[k].y;
       double t_obs_id = transformed_observations[k].id;
       double multi_prob = 1.0;

       for (unsigned int l = 0; l < predicted_landmarks.size(); l++) {
         double pred_landmark_x = predicted_landmarks[l].x;
         double pred_landmark_y = predicted_landmarks[l].y;
         double pred_landmark_id = predicted_landmarks[l].id;

         if (t_obs_id == pred_landmark_id) {
                  double sigma_x = std_landmark[0];
                  double sigma_y = std_landmark[1];
                  multi_prob = (1.0 / (2.0 * M_PI * sigma_x * sigma_y)) * exp(-1.0 * ((pow((t_obs_x - pred_landmark_x), 2) /(2.0 * pow(sigma_x, 2))) + 
                                                                          (pow((t_obs_y - pred_landmark_y), 2)  / (2.0 * pow(sigma_y, 2)))));
                  particles[i].weight *= multi_prob;
         }
       }
     }
    weight_normalizer += particles[i].weight;
   }

   /*Step 5: Normalize the weights of all particles.*/
   for (int i = 0; i < particles.size(); i++) {
      particles[i].weight /= weight_normalizer;
      weights[i] = particles[i].weight;
   }
}
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> resampled_particles;


  // generate random starting index for resampling wheel
  uniform_int_distribution<int> particle_index(0, num_particles - 1);

  int current_index = particle_index(gen);

  double beta = 0.0;
  // get max weight
  double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());
  // spin the resample wheel!
  for (int i = 0; i < particles.size(); i++) {
    // uniform random distribution [0.0, max_weight)
    uniform_real_distribution<double> random_weight(0.0, max_weight_2);
    beta += random_weight(gen);

    while (beta > weights[current_index]) {
        beta -= weights[current_index];
        current_index = (current_index + 1) % num_particles;
    }
      resampled_particles.push_back(particles[current_index]);
  }
  particles = resampled_particles;
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






