function [Y, U] = PlotGenNmr(Vec, fn)
    % first, clips the solution into two (u, v), then form matrix
    % clips a vector in the midpoint of matrix

    global nx ny Lx Ly...
    uBC_L uBC_R uBC_B uBC_T vBC_L vBC_R vBC_T vBC_B
    
    [Mat_u, Mat_v] = VectToMat(Vec);
        % clipping the input vect' into 2, n turning them into mat'
    
    Mat_u(1,:) = uBC_L; % Replacing the BCs within sol'
%     Mat_v(:,1) = vBC_B;
        % This si redundant part thought

    Mat_u_vw = flipud(Mat_u'); % making similar to diagram
%     Mat_v_vw = flipud(Mat_v');
    
    %% Clipping-off u1 values
        % ..........across channel (along x2 line), at channel center
   if mod(nx,2) == 1   % odd num'
       q = (nx+1)/2;    % middle, not CG line
       K1 = Mat_u_vw(:,q);
       K2 = Mat_u_vw(:,q+1); % clipping off the col' across CG
       K_vw = (K1 + K2)/2;
    
    else                % even num'
       q = nx /2;       % middle, not CG line
       K_vw = Mat_u_vw(:, q+1); % clipping off col' across CG
    end
    
    %% Plotter
    u1_y_nm = flipud(K_vw); % reshaping back, for plotting
    u1_y_nm_ex = [uBC_B(q); u1_y_nm; uBC_T(q)];
%     u1_y_nm(ny+1) = uBC_T(q);   % appending the boundary value
    
        %% Grid
    y0 = linspace(-Ly/2, Ly/2, 2*ny +1); % additional 1 for boundaries
%     y_ov = y0(1:2:end-1); % excluding the top boundary
    y_eu = y0(2:2:end); % takes up even values
    y_eu_ex = [-Ly/2, y_eu, Ly/2]; % appended boundary nodes
        
        %% The Plot
    figure(fn)
    plot(y_eu_ex, u1_y_nm_ex, 'LineWidth', 1.5), grid on
    xlabel('Channel Width, $$ x_2(m) $$', 'Interpreter', 'latex')
    ylabel('Velocity, $$ u_1(m/s) $$', 'Interpreter', 'latex')
    title('$$ x_2 $$ vs $$ u_1 $$', 'Interpreter', 'latex')
    
    %% Output
    Y = y_eu_ex; U = u1_y_nm_ex;

end