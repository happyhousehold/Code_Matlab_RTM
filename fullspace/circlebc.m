function [x,y]=circlebc(t,type)
radius = 1;
if type==1
    x = radius*cos(t);   %% primal function
    y = radius*sin(t);
else if type==2
        x = -radius*sin(t);  %% derivative of order one
        y = radius*cos(t);
    else if type==3     %% derivative of order two
            x = -radius*cos(t);
            y = -radius*sin(t);
        else if type==4  %% derivative of order three
                x = radius*sin(t);
                y = -radius*cos(t);
            end
        end
    end

end


%% the rectangle shape %%%%%%%%%%%%%%%%%%
%  radius = 0.5;
% if type==1
%     x = radius*( cos(t).^3 + cos(t) );   %% primal function
%     y = radius*( sin(t).^3 + sin(t) ) + 4;
% else if type==2
%         x = -radius*sin(t).*(  3*cos(t).^2 + 1  );  %% derivative of order one
%         y = radius*cos(t).*(   3*sin(t).^2 + 1 );
%     else if type==3     %% derivative of order two
%             x = -radius*( 3*cos(t).^3 + cos(t) -6*cos(t).*sin(t).^2 );
%             y =  radius*( -(3*sin(t).^3 + sin(t)) + 6*sin(t).*cos(t).^2 ) ;
%         end
%     end
% end
% 
% return

% theta = t;
% p = 5;
% d = 1;
% rho = d*(1+0.2*cos(p*theta));
% drho = -d*p*0.2*sin(p*theta);
% ddrho =-d*p^2*0.2*cos(p*theta);
% if type==1
%     x = rho.*cos(theta)-2;   %% primal function
%     y = rho.*sin(t) + 0;
% else if type==2
%         x = -rho.*sin(t)+drho.*cos(t);  %% derivative of order one
%         y = rho.*cos(t)+sin(t).*drho ;
%     else if type==3     %% derivative of order two
%             x = ddrho.*cos(t) - 2*drho.*sin(t) + rho.*cos(t) ;
%             y = ddrho.*sin(t) + 2*drho.*cos(t) - rho.*sin(t);
%         end
%     end
%     
% end
% return
%% p-leaf
% r = 1;
% if type==1
%     x = r*(cos(t) + 0.2*cos(3*t));
%     y = r*(sin(t) + 0.2*sin(3*t)) +4;
% else if type==2
%         x = -r*(sin(t) + 3*0.2*sin(3*t));
%         y = r*(cos(t) + 3*0.2*cos(3*t)) ;
%     else if type==3
%             x = - r*( cos(t) + 3*3*0.2*cos(3*t));
%             y = - r* ( sin(t) + 3*3*0.2*sin(3*t));
%         end
%     end
% end
% return 