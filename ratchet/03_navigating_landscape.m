%script to plot the path of a float as the ratchet device is moved

%load data
load('./ratchet landscapes/ratchet_sweep_20.mat')

%process data
E_scape = float_min(:,:,size(float_min,3)).';
E_scape = 72*(E_scape-E_scape(:,1));
E_scape = smoothdata(E_scape.','gaussian',24).';

%initialize vectors
turn_num = size(E_scape,1);
theta_num = size(E_scape,2);
beta = zeros(turn_num,1);
winding = zeros(turn_num,1);

%initialize winding at 0 and initial float position parallel to the ratchet channel
min_guess = 1;
w = 0;
for i = 1:turn_num
    
    %perform greedy descent on the energy landscape of the float
    %with a fixed channel geometry.

    %code moves the guess of the minimum downhill in the energy landscape.
    %if the float rotates around in the process, then the winding is updated.

    %when the float can no longer move downhill, we move on to the subsequent
    %channel geometry. the current float position is used to initialize
    %the next round of minimization.
    minimum_searching = true;
    while minimum_searching
        
        if min_guess == 1
            current_index = min_guess;
            ahead_index = min_guess+1;
            behind_index = theta_num;
        elseif min_guess == theta_num
            current_index = min_guess;
            ahead_index = 1;
            behind_index = min_guess-1;
        else
            current_index = min_guess;
            ahead_index = min_guess+1;
            behind_index = min_guess-1;
        end
        
        E0 = E_scape(i,current_index);
        EP = E_scape(i,ahead_index);
        EM = E_scape(i,behind_index);

        if E0-EP > 0
            min_guess = ahead_index;
            if ahead_index == 1
                w = w+1;
            else
            end
        elseif E0-EM > 0
            min_guess = behind_index;
            if behind_index == theta_num
                w = w-1;
            else
            end
        else
            minimum_searching = false;
            beta(i) = min_guess;
            winding(i) = w;
        end
        
    end
end

B = (beta+winding*theta_num)/theta_num*pi;
plot(B)