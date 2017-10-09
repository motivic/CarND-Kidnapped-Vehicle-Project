/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;
	best_particle_idx = -1;

	// Create particles
	default_random_engine gen;	// Set random number generator
	normal_distribution<double> dist_x(x, std[0]);	// Set distribution for x-coordinate
	normal_distribution<double> dist_y(y, std[1]); 	// Set distribution for y-coordinate
	normal_distribution<double> dist_theta(theta, std[2]);	// Set distribution for theta

	Particle particle_tmp;
	for (unsigned int i=0; i<num_particles; ++i) {
		particle_tmp.id = i;
		particle_tmp.x = dist_x(gen);
		particle_tmp.y = dist_y(gen);
		particle_tmp.theta = dist_theta(gen);
		particle_tmp.weight = 1.0/num_particles;
		particles.push_back(particle_tmp);
		weights.push_back(1.0/num_particles);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	for (unsigned int i=0; i<num_particles; ++i) {
		// Handle yaw rate equal to zero case and non-zero case differently
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		double x_new;
		double y_new;
		double theta_new;

		if (fabs(yaw_rate) > 0.0001) {
			x_new = x+velocity/yaw_rate*(sin(theta+yaw_rate*delta_t)-sin(theta));
			y_new = y+velocity/yaw_rate*(cos(theta)-cos(theta+yaw_rate*delta_t));
			theta_new = theta+yaw_rate*delta_t;
		}
		else {
			x_new = x+velocity*delta_t*cos(theta);
			y_new = y+velocity*delta_t*sin(theta);
			theta_new = theta;
		}

		normal_distribution<double> dist_x(x_new, std_pos[0]);
		normal_distribution<double> dist_y(y_new, std_pos[1]);
		normal_distribution<double> dist_theta(theta_new, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
}

std::vector<int> ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, 
												 std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	vector<int> associations;
	associations.resize(predicted.size(), -1);

	// Make sure the predicted nor observations is not empty
	if (predicted.size() == 0 || observations.size() == 0)
		return associations;

	// Find the nearest landmark for each observation
	for (unsigned int i=0; i<predicted.size(); ++i) {
		double x = predicted[i].x;
		double y = predicted[i].y;
		double min_dist = dist(observations[0].x, observations[0].y, x, y);
		int min_idx = 0;

		for (unsigned int j=1; j<observations.size(); ++j) {
			double new_dist = dist(observations[j].x, observations[j].y, x, y);
			if (new_dist < min_dist) {
				min_dist = new_dist;
				min_idx = j;
			}
		}
		associations[i] = min_idx;
	}
	return associations;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	// Map landmarks within range
	vector<LandmarkObs> predicted_landmarks;
	
	// Associate observed landmarks to the possible map landmarks if 
	// a sufficiently accurate localization has been made
	vector<int> associations;

	if (best_particle_idx != -1) {
		double best_x = particles[best_particle_idx].x;
		double best_y = particles[best_particle_idx].y;
		double best_theta = particles[best_particle_idx].theta;
	
		// Find landmarks within sensor range
		LandmarkObs landmark_map;
		for (auto& it : map_landmarks.landmark_list) {
			double distance = dist(best_x, best_y, it.x_f, it.y_f);
			if (distance <= sensor_range) {
				landmark_map.id = it.id_i;	
				landmark_map.x = it.x_f;
				landmark_map.y = it.y_f;
				predicted_landmarks.push_back(landmark_map);
			}
		}

		// Transform observations from VEHICLE's coordinate system into the MAP's coordinate system
		vector<LandmarkObs> observations_map;
		LandmarkObs obs_map;
		for (auto& it : observations) {
			double x_obs = it.x;
			double y_obs = it.y;

			obs_map.x = cos(best_theta)*x_obs-sin(best_theta)*y_obs+best_x;
			obs_map.y = sin(best_theta)*x_obs+cos(best_theta)*y_obs+best_y;

			observations_map.push_back(obs_map);
		}

		associations = dataAssociation(predicted_landmarks, observations_map);
	}

	double total_weight = 0.0;
	for (unsigned int i=0; i<num_particles; ++i) {
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		// Find all map landmarks that are within sensor_range
		if (best_particle_idx == -1) {
			predicted_landmarks.clear();
			LandmarkObs landmark_map;

			for (auto& it : map_landmarks.landmark_list) {
				double distance = dist(x, y, it.x_f, it.y_f);
				if (distance <= sensor_range) {
					landmark_map.id = it.id_i;	
					landmark_map.x = it.x_f;
					landmark_map.y = it.y_f;
					predicted_landmarks.push_back(landmark_map);
				}
			}
		}

		// Skip over if there are no landmarks in range
		if (predicted_landmarks.size() == 0) {
			total_weight += weights[i];
			continue;
		}

		// Skip over if the number of observed landmarks is less than 
		// the number of map landmarks within sensor's range
		if (observations.size() < predicted_landmarks.size()) {
			total_weight += weights[i];
			continue;
		}

		// Transform observations from VEHICLE's coordinate system into the MAP's coordinate system
		vector<LandmarkObs> observations_map(observations.size());
		cout << "Number of observations: " << observations.size() << endl;
		LandmarkObs obs_map;
		for (unsigned int j=0; j<observations.size(); ++j) {
			double x_obs = observations[j].x;
			double y_obs = observations[j].y;

			obs_map.x = cos(theta)*x_obs-sin(theta)*y_obs+x;
			obs_map.y = sin(theta)*x_obs+cos(theta)*y_obs+y;

			observations_map[j] = obs_map;
		}

		// If the map-observation associations haven't been made, make the associations
		if (best_particle_idx == -1)
			associations = dataAssociation(predicted_landmarks, observations_map);

		// Calculate weight for this particle
		// In order to prevent arithmetic underflow, we convert the product to a sum of logarithms and exponentiate
		// Also update the particle's association information
		vector<int> association_ids(predicted_landmarks.size());
		vector<double> sense_x(predicted_landmarks.size());
		vector<double> sense_y(predicted_landmarks.size());

		double log_weight = 0.0;
		double log_denom = log(2.0*M_PI*std_x*std_y);
		for (unsigned int j=0; j<predicted_landmarks.size(); ++j) {
			double x_diff = predicted_landmarks[j].x-observations_map[associations[j]].x;
			double y_diff = predicted_landmarks[j].y-observations_map[associations[j]].y;
			double exponent = -x_diff*x_diff/(2*std_x*std_x)-y_diff*y_diff/(2*std_y*std_y);
			log_weight += exponent-log_denom;

			// Update particle
			association_ids[j] = predicted_landmarks[j].id;
			sense_x[j] = observations_map[associations[j]].x;
			sense_y[j] = observations_map[associations[j]].y;
		}
		particles[i] = SetAssociations(particles[i], association_ids, sense_x, sense_y);

		particles[i].weight = exp(log_weight);
		weights[i] = exp(log_weight);
		total_weight += exp(log_weight);
	}

	// Normalize the weights
	if (total_weight > 0.0001) {
		for (unsigned int i=0; i<num_particles; ++i) {
			particles[i].weight /= total_weight;
			weights[i] /= total_weight;
		}
	}
	else {
		// If total weight is too small, re-initialize the filter
		for (unsigned int i=0; i<num_particles; ++i) {
			particles[i].weight = 1.0/num_particles;
			weights[i] = 1.0/num_particles;
		}
		best_particle_idx = -1;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> resampler (weights.begin(), weights.end());

	vector<Particle> new_particles;
	double total_weight = 0.0;
	for (unsigned int i=0; i<num_particles; ++i) {
		int idx = resampler(gen);
		total_weight += particles[idx].weight;
		new_particles.push_back(particles[idx]);
	}

	double best_weight = 1.0/num_particles; 
	for (unsigned int i=0; i<num_particles; ++i) {
		new_particles[i].weight /= total_weight;
		weights[i] = new_particles[i].weight;
		if (weights[i] > best_weight) {
			best_weight = weights[i];
			best_particle_idx = i;
		}
	}

	cout << "Best particle index: " << best_particle_idx << endl;
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}