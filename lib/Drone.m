classdef Drone < handle
    %% MEMBERS
    properties
        g     % Gravitational Acceleration
        t     % Current Time
        dt    % Delta Time
        tf    % Final Time Simulation

        m     % Mass
        l     % Arm Length
        I     % Inertia Matrix

        x     % State X [X,Y,Z,dX,dY,dZ,phi,theta,psi,p,q,r]
        r     % Position R[X,Y,Z]
        dr    % Velocity V[dx,dy,dz]
        euler % Euler Angels [phi,theta,psi]
        w     % [p,q,r]
        dx    % State Derivative X'
        u     % Signal Control U [T, M1, M2, M3, M4]
        T     % Thrust
        M     % [M1, M2, M3, M4]
    end

    properties

        %% ATTITUDE PARAMETER
        phi_des
        phi_err
        phi_err_prev
        phi_err_sum

        theta_des
        theta_err
        theta_err_prev
        theta_err_sum

        psi_des
        psi_err
        psi_err_prev
        psi_err_sum

        zdot_des
        zdot_err
        zdot_err_prev
        zdot_err_sum

        kP_phi
        kI_phi
        kD_phi

        kP_theta
        kI_theta
        kD_theta

        kP_psi
        kI_psi
        kD_psi

        kP_zdot
        kI_zdot
        kD_zdot

        %% POSISITON PARAMETER
        kP_x
        kD_x

        kP_y
        kD_y

        kP_z
        kD_z

        maxAccelXY
        maxAccelZ
        maxTilt
    end
    %% METHODS
    methods
        %% CONSTRUCTOR
        function obj = Drone(params, initState, InitInputs, gains, simTime)
            obj.g = 9.81;
            obj.t = 0.0;
            obj.dt = 0.01;
            obj.tf = simTime;

            obj.m = params('mass');
            obj.l = params('armLength');
            obj.I = [params('Ixx'),  0,              0; ...
                     0,              params('Iyy'),  0; ...
                     0,              0,              params('Izz')];
            obj.x = initState;
            obj.r = obj.x(1:3);
            obj.dr = obj.x(4:6);
            obj.euler = obj.x(7:9);
            obj.w = obj.x(10:12);
            obj.dx = zeros(12,1);
            obj.u = InitInputs;
            obj.T = obj.u(1);
            obj.M = obj.u(2:4);


            obj.phi_des = 0.0;
            obj.phi_err = 0.0;
            obj.phi_err_prev = 0.0;
            obj.phi_err_sum = 0.0;

            obj.theta_des = 0.0;
            obj.theta_err = 0.0;
            obj.theta_err_prev = 0.0;
            obj.theta_err_sum = 0.0;

            obj.psi_des = 0.0;
            obj.psi_err = 0.0;
            obj.psi_err_prev = 0.0;
            obj.psi_err_sum = 0.0;

            obj.zdot_des = 0.0;
            obj.zdot_err = 0.0;
            obj.zdot_err_prev = 0.0;
            obj.zdot_err_sum = 0.0;

            obj.kP_phi = gains('P_phi');
            obj.kI_phi = gains('I_phi');
            obj.kD_phi = gains('D_phi');

            obj.kP_theta = gains('P_theta');
            obj.kI_theta = gains('I_theta');
            obj.kD_theta = gains('D_theta');

            obj.kP_psi = gains('P_psi');
            obj.kI_psi = gains('I_psi');
            obj.kD_psi = gains('D_psi');

            obj.kP_zdot = gains('P_zdot');
            obj.kI_zdot = gains('I_zdot');
            obj.kD_zdot = gains('D_zdot');



            obj.maxAccelXY = 4.0;
            obj.maxAccelZ = 3.0;
            obj.maxTilt = 20 * pi / 180;
            obj.kP_x = 1.0;
            obj.kD_x = 1.0;
            obj.kP_y = 1.0;
            obj.kD_y = 1.0;
            obj.kP_z = 1.0;
            obj.kD_z = 1.0;
        end

        function state = GetState(obj)
            state = obj.x;
        end

        function obj = EvalEOM(obj)
            bRi = RYY2Rot(obj.euler);
            R = bRi';
            obj.dx(1:3) = obj.dr; %x(4:6)
            obj.dx(4:6) = (1/ obj.m) * ([0; 0; obj.m * obj.g] + R * obj.T * [0; 0; -1]);
            phi = obj.euler(1);
            theta = obj.euler(2);
            obj.dx(7:9) = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
                           0 cos(phi)            -sin(phi);
                           0 sin(phi)*sec(theta) cos(phi)*sec(theta)] * obj.w;
            obj.dx(10:12) = (obj.I) \ (obj.M - cross(obj.w, obj.I * obj.w));
        end

        function obj = UpdateState(obj)
            obj.t = obj.t + obj.dt;
            obj.EvalEOM(); 
            obj.x = obj.x + obj.dx.*obj.dt; % x_k+1 = x_k +  x'k * dt
            obj.r = obj.x(1:3);
            obj.dr = obj.x(4:6);
            obj.euler = obj.x(7:9);
            obj.w = obj.x(10:12);
        end

        %% CONTROLLER
        function obj = AttitudeCtrl(obj, refSig)
         
            obj.phi_des = refSig(1);
            obj.theta_des = refSig(2);
            obj.psi_des = refSig(3);
            obj.zdot_des = refSig(4);

            obj.phi_err = obj.phi_des - obj.euler(1);
            obj.theta_err = obj.theta_des - obj.euler(2);
            obj.psi_err = obj.psi_des - obj.euler(3);
            obj.zdot_err = obj.zdot_des - obj.dr(3);

            %% PHI PID
            obj.u(2) = (obj.kP_phi * obj.phi_err + ...
                        obj.kI_phi * obj.phi_err_sum + ...
                        obj.kD_phi * (obj.phi_err - obj.phi_err_prev) / obj.dt);
                        
            obj.phi_err_sum = obj.phi_err_sum + obj.phi_err;
            obj.phi_err_prev = obj.phi_err;


            %% THETA PID
            obj.u(3) = (obj.kP_theta * obj.theta_err + ...
                        obj.kI_theta * obj.theta_err_sum + ...
                        obj.kD_theta * (obj.theta_err - obj.theta_err_prev) / obj.dt);
            
            obj.theta_err_sum = obj.theta_err_sum + obj.theta_err;
            obj.theta_err_prev = obj.theta_err;

            %% PSI PID
            obj.u(4) = (obj.kP_psi * obj.psi_err + ...
                        obj.kI_psi * obj.psi_err_sum + ...
                        obj.kD_psi * (obj.psi_err - obj.psi_err_prev) / obj.dt);
            
            obj.psi_err_sum = obj.psi_err_sum + obj.psi_err;
            obj.psi_err_prev = obj.psi_err;

            %% Z PID
            obj.u(1) = obj.m * obj.g - ...
                       (obj.kP_zdot * obj.zdot_err + ...
                        obj.kI_zdot * obj.zdot_err_sum + ...
                        obj.kD_zdot * (obj.zdot_err - obj.zdot_err_prev) / obj.dt);
            
            obj.zdot_err_sum = obj.zdot_err_sum + obj.zdot_err;
            obj.zdot_err_prev = obj.zdot_err;


           
            obj.T = obj.u(1);
            obj.M = obj.u(2:4);
            

            
        end
        function obj = PositionCtrl(obj, positionSig)

            %% DESIRED POSITION
            x_des = positionSig(1);
            y_des = positionSig(2);
            z_des = positionSig(3);
        
            %% POSITION ERROR
            ex = x_des - obj.r(1);
            ey = y_des - obj.r(2);
            ez = z_des - obj.r(3);
        
            %% VELOCITY ERROR
            evx = -obj.dr(1);
            evy = -obj.dr(2);
            evz = -obj.dr(3);
        
            %% PD POSITION CONTROLLER
            ax_cmd = obj.kP_x * ex + obj.kD_x * evx;
            ay_cmd = obj.kP_y * ey + obj.kD_y * evy;
            az_cmd = obj.kP_z * ez + obj.kD_z * evz;
        
            %% ACCELERATION SATURATION
            ax_cmd = max(min(ax_cmd, obj.maxAccelXY), -obj.maxAccelXY);
            ay_cmd = max(min(ay_cmd, obj.maxAccelXY), -obj.maxAccelXY);
            az_cmd = max(min(az_cmd, obj.maxAccelZ), -obj.maxAccelZ);
        
            %% DESIRED ACCELERATION
            a_cmd = [ax_cmd;
                     ay_cmd;
                     az_cmd];
        
            %% DESIRED THRUST DIRECTION
            thrust_vec = [0;
                          0;
                          obj.g] - a_cmd;
        
            %% THRUST MAGNITUDE
            T_des = obj.m * norm(thrust_vec);
        
            %% DESIRED BODY Z AXIS IN WORLD FRAME
            b3_des = thrust_vec / norm(thrust_vec);
        
            %% CURRENT YAW
            obj.psi_des = obj.euler(3);
        
            %% DESIRED ROLL / PITCH
            obj.phi_des = atan2( ...
                -b3_des(2), ...
                 b3_des(3));
        
            obj.theta_des = atan2( ...
                 b3_des(1), ...
                 sqrt(b3_des(2)^2 + b3_des(3)^2));
        
            %% ANGLE SATURATION
            obj.phi_des = max(min(obj.phi_des, obj.maxTilt), -obj.maxTilt);
            obj.theta_des = max(min(obj.theta_des, obj.maxTilt), -obj.maxTilt);
        
            %% SEND COMMAND TO ATTITUDE CONTROLLER
            obj.phi_des = obj.phi_des;
            obj.theta_des = obj.theta_des;
            obj.psi_des = obj.psi_des;
        
            %% Z VELOCITY COMMAND
            obj.zdot_des = -obj.dr(3);
        
            %% THRUST
            obj.T = T_des;
        
        end
        
    end
end