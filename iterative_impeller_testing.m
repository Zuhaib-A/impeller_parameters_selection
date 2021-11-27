clear
clc

format longg
flow_rate=0.00125;
pa_gradient=80000;
rpm=input("Enter the RPM (2400, 2500, 2600, 2700, 2800) ");
r1=0.017;
r2=0.05;
spiral_angle=5.3;
density=1000;
g=9.81;
leakage_percentage=30;
pressure_loss_percentage=50;
diffuser_area=0.000278;

adjusted_flow_rate=flow_rate*((100+leakage_percentage)/100);
adjusted_pa_gradient=pa_gradient*((100+pressure_loss_percentage)/100);
operating_speed=rpm*(2*pi/60);
net_head=adjusted_pa_gradient/(g*density);

inner_blade_thick=[0.001:0.0001:0.0025];
outer_blade_thick=[0.001:0.0001:0.0025];

b1=[0.008:0.0001:0.019];
b2=[0.005:0.0001:0.0079];
n=input("Enter the number of blades: ");

results_bf1=[];

guessed_bf1=[0.1:0.01:0.95];
guessed_bf2=[0.1:0.01:0.95];

index=0;
for tk=inner_blade_thick
for c=guessed_bf1
for i=b1
    index=index+1; 
    for x=n
        for t=1:1000
            current_beta1_from_guess=(atan((adjusted_flow_rate/(2*pi*i*c*operating_speed*r1^2)))*(180/pi));
            current_calculated_bf1=(1-((x*tk*i)/(2*pi*r1*i*sin(current_beta1_from_guess*(pi/180)))));
            c=current_calculated_bf1;
        end
        checking_c=(1-((x*tk*i)/(2*pi*r1*i*sin(current_beta1_from_guess*(pi/180)))));
       
        if abs(c-checking_c)<=0.01
            
            if c < 1
                if current_beta1_from_guess>=0
                    results_bf1(1,index)=r1;
                    results_bf1(2,index)=current_calculated_bf1;
                    results_bf1(3,index)=i;
                    results_bf1(4,index)=x;
                    results_bf1(5,index)=tk;
                    results_bf1(6,index)=current_beta1_from_guess;      
                    index=index+1;
                else
                end
            else
            end
        else
        end
    end
        index=index-1;
end
end
end

results_bf2=[];
outflow_disalign=[];
real_magnitude_out=[];

index=0;
for tk=outer_blade_thick
for c=guessed_bf2
for i=b2
    index=index+1;
    for x=n
        for t=1:1000
            current_beta2_from_guess=(atan(((adjusted_flow_rate*operating_speed*r2)/(2*pi*r2*i*c*((operating_speed^2*r2^2)-g*net_head))))*(180/pi));
            current_calculated_bf2=(1-((x*tk*i)/(2*pi*r2*i*sin(current_beta2_from_guess*(pi/180)))));
            c=current_calculated_bf2;
        end
        
        checking_c=(1-((x*tk*i)/(2*pi*r2*i*sin(current_beta2_from_guess*(pi/180)))));
        
        if abs(c-checking_c)<=0.05
           
            outlet_normal_vel=(adjusted_flow_rate/(2*pi*r2*i*c));
            outlet_tangent_vel=(g*net_head)/(operating_speed*r2);
            alpha_s=(atan(outlet_normal_vel/outlet_tangent_vel)*(180/pi));
            slip_factor=(1-(((sqrt(sin(current_beta2_from_guess*(pi/180)))))/(x^(0.7))));
            real_outlet_tangent_vel=outlet_tangent_vel*slip_factor;
            real_alpha_s=(atan(outlet_normal_vel/real_outlet_tangent_vel)*(180/pi));
            outflow_disalign=[outflow_disalign, spiral_angle-real_alpha_s];
            real_magnitude_out=[real_magnitude_out, sqrt((outlet_normal_vel)^2+(real_outlet_tangent_vel)^2)];

            if c < 1
                if current_beta2_from_guess > 0
                    results_bf2(1,index)=r2;
                    results_bf2(2,index)=current_calculated_bf2;    
                    results_bf2(3,index)=i;
                    results_bf2(4,index)=x;
                    results_bf2(5,index)=tk;
                    results_bf2(6,index)=current_beta2_from_guess;

                    results_bf2(7,index)=outflow_disalign(index);
                    results_bf2(8,index)=real_magnitude_out(index)-(adjusted_flow_rate/diffuser_area);
                    index=index+1;
                else
                end
        else
        end
        else
        end
    end
        index=index-1;
end
end
end

results_bf1 = unique(results_bf1.', 'rows', 'stable').';
results_bf2 = unique(results_bf2.', 'rows', 'stable').';

[value2, index_bf_outflow]=(mink(abs(((results_bf2(7,:))'))',10));
[value3, index_velocity_diff]=(mink(abs(((results_bf2(8,:))'))',10));
[value4, index_blockage1]=(maxk(abs(((results_bf1(2,:))'))',10));
[value5, index_blockage2]=(maxk(abs(((results_bf2(2,:))'))',10));

best_results_outflow_bf2=[];
index=0;
for x = index_bf_outflow
    index=index+1;
    best_results_outflow_bf2(1,index)=results_bf2(1,x);
    best_results_outflow_bf2(2,index)=results_bf2(2,x);
    best_results_outflow_bf2(3,index)=results_bf2(3,x);
    best_results_outflow_bf2(4,index)=results_bf2(4,x);
    best_results_outflow_bf2(5,index)=results_bf2(5,x);
    best_results_outflow_bf2(6,index)=results_bf2(6,x);
    best_results_outflow_bf2(7,index)=results_bf2(7,x);
    best_results_outflow_bf2(8,index)=results_bf2(8,x);
end
best_results_velocity_bf2=[];
index=0;
for x = index_velocity_diff
    index=index+1;
    best_results_velocity_bf2(1,index)=results_bf2(1,x);
    best_results_velocity_bf2(2,index)=results_bf2(2,x);
    best_results_velocity_bf2(3,index)=results_bf2(3,x);
    best_results_velocity_bf2(4,index)=results_bf2(4,x);
    best_results_velocity_bf2(5,index)=results_bf2(5,x);
    best_results_velocity_bf2(6,index)=results_bf2(6,x);
    best_results_velocity_bf2(7,index)=results_bf2(7,x);
    best_results_velocity_bf2(8,index)=results_bf2(8,x);
end
best_results_blockage1=[];
index=0;
for x = index_blockage1
    index=index+1;
    best_results_blockage1(1,index)=results_bf1(1,x);
    best_results_blockage1(2,index)=results_bf1(2,x);
    best_results_blockage1(3,index)=results_bf1(3,x);
    best_results_blockage1(4,index)=results_bf1(4,x);
    best_results_blockage1(5,index)=results_bf1(5,x);
    best_results_blockage1(6,index)=results_bf1(6,x);
end
best_results_blockage2=[];
index=0;
for x = index_blockage2
    index=index+1;
    best_results_blockage2(1,index)=results_bf2(1,x);
    best_results_blockage2(2,index)=results_bf2(2,x);
    best_results_blockage2(3,index)=results_bf2(3,x);
    best_results_blockage2(4,index)=results_bf2(4,x);
    best_results_blockage2(5,index)=results_bf2(5,x);
    best_results_blockage2(6,index)=results_bf2(6,x);
    best_results_blockage2(7,index)=results_bf2(7,x);
    best_results_blockage2(8,index)=results_bf2(8,x);
end

display('Assuming a leakage percentage of '+ string(leakage_percentage) + '% and a pressure loss of ' + string(pressure_loss_percentage)+ '%, the best options are as follows:')

rows={'r1','inner blockage factor','b1', 'number of blades','inner blade thickness', 'beta 1'};
cols=["1","2","3", "4", "5", "6", "7", "8", "9", "10"];
best_results_blockage1_tb=array2table(best_results_blockage1,'RowNames',rows,'VariableNames',cols)

rows={'r2','outer blockage factor','b2', 'number of blades','outer blade thickness', 'beta 2', 'outflow disalign', 'outlet velocity difference'};
cols=["1","2","3", "4", "5", "6", "7", "8", "9", "10"];
best_results_outflow_bf2_tb=array2table(best_results_outflow_bf2,'RowNames',rows,'VariableNames',cols)

rows={'r2','outer blockage factor','b2', 'number of blades','outer blade thickness', 'beta 2', 'outflow disalign', 'outlet velocity difference'};
cols=["1","2","3", "4", "5", "6", "7", "8", "9", "10"];
best_results_velocity_bf2_tb=array2table(best_results_velocity_bf2,'RowNames',rows,'VariableNames',cols)

rows={'r2','outer blockage factor','b2', 'number of blades','outer blade thickness', 'beta 2', 'outflow disalign', 'outlet velocity difference'};
cols=["1","2","3", "4", "5", "6", "7", "8", "9", "10"];
best_results_blockage2_tb=array2table(best_results_blockage2,'RowNames',rows,'VariableNames',cols)