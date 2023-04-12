#ifndef GET_RID_OF_TRANSITION
#define GET_RID_OF_TRANSITION

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>


class GetRidofTransition
{
        public:
                GetRidofTransition(double timestep, int delay_passages_n,std::vector<double>& traj, int N_position, int N_velocity)
                        : timestep_(timestep),
                          traj_(traj),
                          delay_passages_N(delay_passages_n),
                          N_position_bins(N_position),
                          N_velocity_bins(N_velocity)
        {}
                void Locate_transition();
                //std::vector<std::vector<double>> Without_transition_traj;

                std::vector<std::vector<long long int>> heatmap;
                // Build Velocity-Position histogram
                void Set_bins(); 
                void Save_data(std::string filename);




        private:
                double timestep_;
                std::vector<double> traj_;
                std::vector<double> velocity;
                double initial_state;
                double Last_state;
                std::vector<double> time_passage_up;
                std::vector<double> time_passage_down;
                std::vector<std::vector<int>> downcrossing_index;
                std::vector<std::vector<int>> upcrossing_index;
                double delay_passages_N;

                // State boundary
                static constexpr double upboundary  = 1.0;
                static constexpr double downboundary = -1.0;
                
                // Position bin
                int N_position_bins;
                std::vector<double> position_bins;
                // Velocity bin
                int N_velocity_bins;
                std::vector<double> Velocity_bins;
                // Heatmap for velocity bins.
                // These two vectors save the bin index for velocity and position
                std::vector<int> Rescale_position;
                std::vector<int> Rescale_velocity;
                std::vector<long long int> heattemp;
};

void GetRidofTransition::Locate_transition()
{
        /* Locate the transition 
         *      We record each crossing to boundaries
         *      1. Recall the initial state
         */
        /*
         * The idea for this part:
         * Since we got each index of the
         * ----|-------------------|----upcrossing_index
         *     .\                  |
         *     . \     /\    /\    |
         * -------\---/--\--/--\--/------downcrossing_index
         *     .  .\ /.  .\/.  .\/..
         * corresponding index to position
         *     .  .   .  .  .  .  ..
         *    up[0]   .  .  .  .  .up[1]
         *        .   .  .  .  .  .
         *   down[0]  1  2  3  4  5
         * The passage time in the lower well will be each index minus the next first element in
         * upcrossing_index.
         */
 
        // The absolute value indicate if there is crossing
        double state = 0;
        double cross_state = 0;
        int jump =0;
        int beginning_state;
        int beginning_point;
        std::cout<<"Find initial State and crossing state"<<std::endl;
        for(unsigned int i; i<traj_.size(); i++){
                if(traj_[i] < downboundary)
                {
                        state = downboundary;
                        beginning_state = 0;
                        cross_state = 0;
                        break;
                }
                if(traj_[i] > upboundary)
                {
                        state = upboundary;
                        beginning_state = 1;
                        cross_state = 1;
                        break;
                }
        }
        //For temperary saving vector, then push_back.
        std::vector<int> uptemp;
        std::cout<<"initial state is "<<beginning_state<<std::endl;
        std::vector<int> downtemp;
 
        // Indentify the passage and transitions
        for(unsigned int i=beginning_point; i<traj_.size(); i++)
        {
                // System starts from downboundary
                if(state == downboundary)
                {
                        //If the state flips, should clear the temp vector and record the first crossing
                        //events
                        if(traj_[i]>=upboundary){
                                //jump to another satate
                                state = upboundary;
	                        cross_state = 1;
                                jump += 1;
                                uptemp.push_back(i);
                                // For passage time test
                                downcrossing_index.push_back(downtemp);
                                downtemp.clear();
                        }
                        else
                        {
                                if(traj_[i]<downboundary & cross_state == 1)
                                {
                                        cross_state = 0;
                                        downtemp.push_back(i);
                                }
                                if(traj_[i]>downboundary & cross_state == 0 )
                                {
                                        cross_state = 1;
                                        downtemp.push_back(i);
                                }
                        }

                }
                else
                {
                        if(traj_[i]<=downboundary){
                                state = downboundary;
                                cross_state = 0;
                                jump += 1;
                                downtemp.push_back(i);

                                upcrossing_index.push_back(uptemp);
                                uptemp.clear();
                        }
                        else
                        {
                                if(traj_[i]<upboundary & cross_state == 1)
                                {
                                        cross_state = 0;
                                        uptemp.push_back(i);
                                }
                                if(traj_[i]>upboundary & cross_state == 0 )
                                {
                                        cross_state = 1;
                                        uptemp.push_back(i);
                                }
                        }
                }



        }
        


        std::cout<<"Testing vector"<<std::endl;
        for(size_t i=0; i<downcrossing_index.size(); i++){
                for(auto & p : downcrossing_index[i]){
                        if(p<0)
                                std::cout<<"Things wrong here"<<p<<std::endl;
                }
        }
        std::cout<<"State index recording end!"<<std::endl;
        // Delete  the transition path and a range of relax time (Let the particles relax).
        // ****************
        // Input is delay_passages_N
        // ****************_
        // For the correlation decay time, we use the passages time
        // likewise
        // 5 passages before the transition and 5 passages after the transition.
        // Which means this relax time remains unclear, more depends on the instinction
        // or experience.
 
        // Relaxing passages time;
        int Relax_ = delay_passages_N;
        if(jump<2){
                std::cout<<"Transition occurs only once, no statistic inside"<<std::endl;
                exit(1);
        }
        std::cout<<"This is the first element:"<<downcrossing_index[18][0]<<std::endl;
        std::cout<<"This is the last element"<<downcrossing_index[18].back()<<std::endl;
        /*
        std::cout<<"This is the wrong index"<<downcrossing_index[18][Relax_]+66670500<<std::endl;
        std::cout<<"This is the wrong index"<<Rescale_position[downcrossing_index[18][Relax_]+66670500]<<std::endl;
        std::cout<<"This is the wrong index"<<Rescale_velocity[downcrossing_index[18][Relax_]+66670500]<<std::endl;
        std::cout<<"This is the wrong index"<<downcrossing_index[18][Relax_+66670544]<<std::endl;
        std::cout<<"This is the wrong index"<<Rescale_position[downcrossing_index[18][Relax_]+66670544]<<std::endl;
        std::cout<<"This is the wrong index"<<Rescale_velocity[downcrossing_index[18][Relax_]+66670544]<<std::endl;
        */
        int size_of_up;
        int size_of_down;
        std::cout<<"Jump time:"<<jump/2<<std::endl;
        int count_index = 0;
        for(int i=0; i<jump/2 -1; i++)
        {
                count_index = 0;
                if(upcrossing_index[i].size()>2*Relax_+2)
                {
                        size_of_up = upcrossing_index[i].size();
                        std::cout<<"Up Start from index: "<<upcrossing_index[i][Relax_]<<std::endl;

                        for(size_t j=upcrossing_index[i][Relax_]; j<upcrossing_index[i][size_of_up-Relax_];j++)
                        {
                                heatmap[Rescale_position[j]][Rescale_velocity[j]] +=1;
                        }
                        std::cout<<"Up End index: "<<upcrossing_index[i][size_of_up-Relax_]<<std::endl;
                }

                if(downcrossing_index[i].size()>2*Relax_+2)
                {
                        size_of_down = downcrossing_index[i].size();
                        for(size_t j=downcrossing_index[i][Relax_]; j<downcrossing_index[i][size_of_down-Relax_];j++)
                        {
                                if(count_index>66670000 | count_index%1000==0){
                                        //std::cout<<"Counts: "<<count_index<<std::endl;
                                }
                                count_index += 1;
                                if(Rescale_position[j]>98 | Rescale_velocity[j]>98)
                                {
                                        std::cout<<"Down End from index: "<<downcrossing_index[i][size_of_down]<<std::endl;
                                }
                                if(Rescale_position[j]<0 | Rescale_velocity[j]<0)
                                {
                                        std::cout<<"Index less than 0"<<std::endl;
                                }
                                heatmap[Rescale_position[j]][Rescale_velocity[j]] +=1;
                        }
                }
                std::cout<<"Circle crashed: "<<i<<std::endl;
        }
        std::cout<<"Transfer index to time"<<std::endl;
        // For calculate the mean all first passage time;
        std::cout<<"Boundary is "<<upboundary<<" "<<downboundary;
}




void GetRidofTransition::Set_bins()
{
        // Set position bins
        double max_traj = *std::max_element(traj_.begin(), traj_.end());
        double min_traj = *std::min_element(traj_.begin(), traj_.end());

        std::cout<<"X Max: "<<max_traj<<std::endl;
        std::cout<<"X Min:"<<min_traj<<std::endl;
        int X_min = std::floor(min_traj);
        int X_max = std::ceil(max_traj);
        if (-min_traj>max_traj)
        {
                X_max = -X_min;
        }
        else
        {
                X_min = -X_max;}


        // Set the symmetry boundary
        double size_of_position_bins = (X_max-X_min)/(double)N_position_bins;
        // Because the boundary, need to add 2 additional data points
        for(int i=0; i<N_position_bins+1; i++)
        {
                position_bins.push_back(X_min + i*size_of_position_bins);
        }
        double Rescale_X = N_velocity_bins/(double)(X_max-X_min);
        X_max = std::ceil(X_max*Rescale_X);
        X_min = std::floor(X_min*Rescale_X);

        // Initialized the velocity
        std::cout<<"Position bins Set"<<std::endl;
        // Get each point velocity (middle point)
        for(unsigned int i=0; i<traj_.size()-2; i++)
        {
                // Record the velocity
                velocity.push_back((traj_[i+2]-traj_[i])/timestep_/2);
                // Initialize Poisition
                // Record the corresponding position (index, not the exact value)
                Rescale_position.push_back(static_cast<int>(std::floor(traj_[i+1]*Rescale_X))-X_min); 
        }

        /* Build velocity Bins
         * 1. Find the maximum and minimum velocity;
         * 2. Rescale the velocity;
        */
        std::cout<<"Velocity length"<<velocity.size()<<std::endl;
        double max_velocity = *std::max_element(velocity.begin(), velocity.end());
        double min_velocity = *std::min_element(velocity.begin(), velocity.end());


        std::cout<<"V Max: "<<max_velocity<<std::endl;
        std::cout<<" V Min:"<<min_velocity<<std::endl;
        int V_min = std::floor(min_velocity);
        int V_max = std::ceil(max_velocity);
        if (-min_velocity>max_velocity)
        {
                V_max = -V_min;
        }
        else
        {
                V_min = -V_max;}

        // Set the symmetry boundary
        double size_of_velocity_bins = 2*V_max/(double)N_velocity_bins;
        double Rescale_V = N_velocity_bins/(double)(V_max-V_min);
        // Because the boundary, need to add 2 additional data points
        for(int i=0; i<N_velocity_bins+1; i++)
        {
                Velocity_bins.push_back(V_min + i*size_of_velocity_bins);
        }

        V_max = std::ceil(V_max*Rescale_V);
        V_min = std::floor(V_min*Rescale_V);
        // Initialized the velocity
        for(unsigned int i=0; i<velocity.size(); i++)
        {
                Rescale_velocity.push_back(static_cast<int>(std::floor(velocity[i]*Rescale_V)-V_min));
        }
        // Initialize the Heatmap;
        std::cout<<"Initializing heatmap"<<std::endl;
        for(int j=0; j<N_velocity_bins; j++){
                heattemp.push_back(0);
        }
        std::cout<<"Bins number of Velocity: "<<heattemp.size();
        for(int i=0; i<N_position_bins; i++)
        {
               heatmap.push_back(heattemp);
        }
        heattemp.clear();
        std::cout<<"Initializing ends"<<std::endl;
}

void GetRidofTransition::Save_data(std::string filename)
{
        // Save Velocity_bins
        std::ofstream V_BINS;
        V_BINS.open("./V_bins.txt");
        for(int v=0; v<Velocity_bins.size(); v++)
        {
                V_BINS<<Velocity_bins[v]<<std::endl;
        }
        V_BINS.close();

 
        // Save Position bins
        std::ofstream X_BINS;
        X_BINS.open("./X_bins.txt");
        for(int x=0; x<position_bins.size(); x++)
        {
                X_BINS<<position_bins[x]<<std::endl;
        }
        X_BINS.close();


        // Save heatmap
        std::ofstream HEATMAP_OUT;
        HEATMAP_OUT.open(filename);
        for(int x=0; x<N_position_bins; x++)
        {
                for(int v=0; v<N_velocity_bins; v++)
                {
                        HEATMAP_OUT<<heatmap[x][v]<<" ";
                }
                HEATMAP_OUT<<std::endl;
        }
        HEATMAP_OUT.close();
}
#endif // GET_RID_OF_TRANSITION


