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

std::map<double,double> sin_theta_map;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;

	random_device rd;
	mt19937 gen(rd());
	// default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);

	particles.resize(num_particles);
	weights.resize(num_particles) ;

	for (int i = 0; i < num_particles; i++) {

		particles[i].id = i;
		particles[i].x = dist_x(gen);  // Add random Gaussian noise to each particle.
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_psi(gen);
		particles[i].weight = 1.0;
		// cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_psi << endl;
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	random_device rd;
	mt19937 gen(rd());

	for (int i = 0; i < num_particles; i++) {

		const double theta = particles[i].theta;
		const double sin_theta = sin(theta);
		const double cos_theta = cos(theta);
		// if yaw_rate is near 0 we are going straight.
		if (fabs(yaw_rate) > 0.00001) {		
			// So we are turning
			const double velocity_yaw_rate = velocity / yaw_rate;
			const double yaw_rate_delta_t = yaw_rate * delta_t;
			particles[i].x += velocity_yaw_rate * (sin(theta + yaw_rate_delta_t) - sin_theta);
			particles[i].y += velocity_yaw_rate * (cos_theta - cos(theta + yaw_rate_delta_t));
			particles[i].theta += yaw_rate * delta_t;
			// Calculate noise on theta
			normal_distribution<double> dist_psi(particles[i].theta, std_pos[2]);
			particles[i].theta = dist_psi(gen);
		}

		else {		
			// we are not turning. Don´t calculate theta.
			const double velocity_delta_t = velocity * delta_t;
			particles[i].x += velocity_delta_t * cos_theta;
			particles[i].y += velocity_delta_t * sin_theta;
		}
		// add random Gaussian noise to x and y. Common to straight and turning.
		normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
		normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
		particles[i].x = dist_x(gen);  
		particles[i].y = dist_y(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

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

	// The firs, some precalculation
	const double standard_x = std_landmark[0];
	const double standard_y = std_landmark[1];
	const double standard_x_squared_2 = (standard_x * standard_x);
	const double standard_y_squared_2 = (standard_y * standard_y);
	const long double w_common_part = 1 / (2 * M_PI * standard_x * standard_y);	

	long double particle_weight = 1.0;
	double threshold = 50.0;
	// A variable for the id of the landmark that is nearest neighbor
	long n_n_id = -1;
	// Declaration of variables to be used on the weight phase to avoid redeclaration on each loop step.
	long double x;
	long double u_x;
	long double y;
	long double u_y;
	long double x_ux_squared;
	long double y_uy_squared;
	long double landmark_part;
	double difference;

	for (int i = 0; i < num_particles; i++) {
		// Precalculation of particle data.
		const double particle_theta = particles[i].theta;
		const double cos_theta = cos(particle_theta);
		const double sin_theta = sin(particle_theta);
		particle_weight = 1.0;
		n_n_id=-1;
		for (unsigned int j = 0; j < observations.size(); j++) {
			// For each observation calculate the coordinates in map reference.
			double tr_obs_x = particles[i].x + (observations[j].x * cos_theta) -
				(observations[j].y * sin_theta);            
			double tr_obs_y = particles[i].y + (observations[j].x * sin_theta) +
				(observations[j].y * cos_theta);
			// Initialize treshold	
			threshold = 50.0;
			// originally, observations has no id, we have to give one id to the best landmark.
			for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++) {

				// Find the predicted measurement that is closest to each observed measurement
				difference = dist(tr_obs_x, tr_obs_y,
					map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
				// check if threshold is less than previous
				if (difference < threshold) {
					// if it is, then we have the nearest neighbour ID
					// map data is not zero indexed
					n_n_id = map_landmarks.landmark_list[k].id_i - 1;
					// set the diff as a new threshold to compare with
					threshold = difference;
				}
			}
			// If there is a good landmark as neighbor. Calculate the weight.
			if (n_n_id >= 0) {

				x = tr_obs_x;
				u_x = map_landmarks.landmark_list[n_n_id].x_f;
				y = tr_obs_y;
				u_y = map_landmarks.landmark_list[n_n_id].y_f;
				x_ux_squared = (x - u_x) * (x - u_x);
				y_uy_squared = (y - u_y) * (y - u_y);
				landmark_part = x_ux_squared / (standard_x_squared_2)+
					y_uy_squared / (standard_y_squared_2);

				//cout << "transformed_observations\t" << x << "," << y << endl ;
				//cout << "\t reasonable_landmarks\t" << u_x << "," << u_y << endl ;
				
				particle_weight *= w_common_part * exp((-1 / 2.) * (landmark_part));

				//cout << particle_weight << endl;
			}
		}
		weights[i]      = particle_weight ;
		particles[i].weight = particle_weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	random_device rd;
	mt19937 gen(rd());
	discrete_distribution<> d(weights.begin(), weights.end());

	vector<Particle> resampled_particles(num_particles);

	for (int i = 0; i < num_particles; i++) {
		resampled_particles[i] = particles[d(gen)];
	}

	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
