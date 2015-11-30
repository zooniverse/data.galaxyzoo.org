
// begin Parameters
function Parameters()
{
    this.potential_type = 1;
    this.ndim = 3;

    this.mass1 = 0;
    this.epsilon1 = 0;
    this.rin1 = 0;
    this.rout1 = 0;
    this.rscale1 =  new Array(3);
    this.theta1=0;
    this.phi1=0;
    this.heat1=0;
    this.opt1=0;

    this.mass2=0;
    this.epsilon2=0;
    this.rin2=0;
    this.rout2=0;
    this.rscale2 = new Array(3);
    this.theta2=0;
    this.phi2=0;
    this.heat2=0;
    this.opt2=0;

    this.eps1=0;
    this.eps2=0;

    this.tcurrent=0;


//  real (kind=8), dimension(:,:), allocatable :: x0, xout
    this.n = 0;
    this.n1 = 0;
    this.n2 = 0;

    this.time = 0;
    this.tstart=0;
    this.tend=0;
    this.tIsSet=false;
    this.inclination_degree=0;
    this.omega_degree=0;
    this.rmin=0;
    this.velocity_factor=0;
    this.mratio=0;
    this.secondary_size=0;

    this.sec_vec = new Array(6);
    this.use_sec_vec = false;

    var hbase = 0.001;
    this.h = 0.001;

    this.nstep = 0;
    this.nout = 0;
    this.iout = 0;
    this.unit = 0;
    this.istep = 0;

    this.fname = "";
    this.iostat = 0;
    this.header_on=false;

    // info about time of closest approach
    this.tmin = 0;
    this.vrmin = 0;
    this.beta = 0;

    this.standardGalaxy1 = function()
    {
        // disk profile #1
        this.mass1 = 1.0;
        this.epsilon1  = 0.3;
        this.rin1 = 0.05;
        this.rout1 = 1.0;
        this.rscale1[0] = 3.0;
        this.rscale1[1] = 3.0;
        this.rscale1[2] = 3.0;
        this.theta1 = 0.0;
        this.phi1 = 0.0;
        this.opt1 = 1;
        this.heat1 = 0.0;
        this.eps1 = this.epsilon1*this.epsilon1;
    }

    this.standardGalaxy2 = function()
    {
        // disk profile #2
        this.mass2 = 1.0;
        this.epsilon2  = 0.3;
        this.rin2 = 0.05;
        this.rout2 = 1.0;
        this.rscale2[0] = 3;
        this.rscale2[1] = 3;
        this.rscale2[2] = 3;
        this.theta2 = 0.0;
        this.phi2 = 0.0;
        this.opt2 = 1;
        this.heat2 = 0.0;
        this.eps2 = this.epsilon2*this.epsilon2;
    }

    this.testCollision = function()
    {
        // collision parameters
        this.inclination_degree = 90.0;
        this.omega_degree = 0.0;
        this.rmin = 1.0;
        this.velocity_factor = 1.0;
        this.time = -3.0;

        // time step
        this.h = 0.001;//hbase;
        this.nout = 5;
        this.nstep = 500;

        // particle numbers
        this.n1 = 1000;
        this.n2 = 1000;
        this.n = this.n1 + this.n2;
    }

    this.customCollision = function()
    {
        this.phi1 = 5.0;
        this.theta1 = 5.0;
        this.rscale1[0] = 1.0;
        this.rscale1[1] = 1.0;
        this.rscale1[2] = 1.0;
        this.rout1 = 1.0;
        this.mass1 = 1.0;
        this.epsilon1 = 0.3;
        this.eps1 = this.epsilon1*this.epsilon1;
        this.n1 = 1000;
        this.heat1 = 0.0;
        this.opt1 = 1;

        this.phi2 = 0.0;
        this.theta2 = 0.0;
        this.rscale2[0] = 0.3;
        this.rscale2[1] = 0.3;
        this.rscale2[2] = 0.3;
        this.rout2 = 0.5;
        this.mass2 = 0.5;
        this.epsilon2 = 0.3;
        this.eps2 = this.epsilon2*this.epsilon2;
        this.n2 = 500;
        this.heat2 = 0.0;
        this.opt2 = 1;

        this.inclination_degree = 20.0;
        this.omega_degree = 0.0;
        this.rmin = 0.9;
        this.velocity_factor = 0.9;

        this.h = 0.001;//hbase;
        this.time = -5;

        this.n = this.n1 + this.n2;
    }

    this.defaultParameters = function()
    {
        this.standardGalaxy1();
        this.standardGalaxy2();
        this.testCollision();

        this.customCollision();
    }
}
// end Parameters



/**
 * This ForceModel represents a softened point-mass.
 *
 */
function SPMModel()
{
    // temporary storage arrays dependent on n
    this.r22;
    this.r21;
    this.r1;
    this.r2;
    this.rn;
    this.a1;
    this.a2;
    this.a3;

    // temporary storage independent of n
    this.xn = new Array(6);
    this.xp = 0.0;
    this.yp = 0.0;
    this.zp = 0.0;

    // force parameters
    this.params = null;
    this.m1;
    this.m2;
    this.m3;

    this.eps1;
    this.eps2;

    /**
     * Sets the simulation parameters to this model.
     */
    this.setParameters = function(p)
    {
        this.params = p;

        if(this.params != null)
        {
            this.m1 = this.params.mass1;
            this.m2 = this.params.mass2;
            this.m3 = this.params.mass2;
            this.eps1 = this.params.epsilon1 * this.params.epsilon1;
            this.eps2 = this.params.epsilon2 * this.params.epsilon2;
        }
    }


    /** 
     * Initialize temporary storage.
     */
    this.initVars = function(n)
    {
        this.r22 = new Array(n);
        this.r21 = new Array(n);
        this.r1 = new Array(n);
        this.r2 = new Array(n);
        this.rn = new Array(n);
        this.a1 = new Array(n);
        this.a2 = new Array(n);
        this.a3 = new Array(n);
        for(var i=0; i<n; i++)
        {
            this.r22[i]=0;
            this.r21[i]=0;
            this.r1[i]=0;
            this.r2[i]=0;
            this.rn[i]=0;
            this.a1[i]=0;
            this.a2[i]=0;
            this.a3[i]=0;
        }
    }

    /** 
     * Cleanup temporary storage.
     */
    this.deallocateVars = function()
    {
        this.r22 = null;
        this.r21 = null;
        this.r1 = null;
        this.r2 = null;
        this.rn = null;
        this.a1 = null;
        this.a2 = null;
        this.a3 = null;
    }

    /**
     * Compute the circular velocity for a particle at a distance r from the specified mass.
     * The rout scale of the disk and softening length, eps, are provided.
     */
    this.circularVelocity = function(mass, r, rout, eps)
    {
        var ftotal = mass / ( r*r + eps );
        var v = Math.sqrt(ftotal*r);

        return v;
    }

    /**
     * For the given particle positions and velocities, calculate the accelerations.
     */
    this.diffeq = function(x, f)
    {
        var tmp1 = 0;
        var tmp2 = 0;
        var tmp3 = 0;
        var tmp4 = 0;
        var tmp5 = 0;

        var n = x.length;
        var i;

        for(i=0;i<6;i++)
        {
            this.xn[i] = x[n-1][i];
        }

        // make temps to handle perturber galaxy
        this.xp = this.xn[0];
        this.yp = this.xn[1];
        this.zp = this.xn[2];

        tmp4 = this.xp*this.xp+this.yp*this.yp+this.zp*this.zp;   // r2n
        tmp5 = Math.sqrt(tmp4);     // rn
        tmp4 = -this.m3/(tmp4+this.eps2);

        for(i=0;i<n;i++)
        {
            tmp1 = (x[i][0]-this.xp);
            tmp2 = (x[i][1]-this.yp);
            tmp3 = (x[i][2]-this.zp);

            this.r22[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;

            tmp1 = x[i][0];
            tmp2 = x[i][1];
            tmp3 = x[i][2];

            this.r21[i] = tmp1*tmp1+tmp2*tmp2 + tmp3*tmp3;

            this.r2[i]  = Math.sqrt(this.r22[i]);
            this.r1[i]  = Math.sqrt(this.r21[i]);
            this.rn[i]  = tmp5;

            this.a1[i] = -this.m1 / (this.r21[i] + this.eps1);
            this.a2[i] = -this.m2 / (this.r22[i] + this.eps2);
            this.a3[i] = tmp4;
        }

        // this is a correction to prevent NaN errors in the vectorized
        // function evalution at the location of the second mass
        this.r2[n-1] = 1.0;

        // calculate the RHS of the diffeq

        // f(:,1) = x(:,4)
        // f(:,2) = x(:,5)
        // f(:,3) = x(:,6)

        // f(:,4) = a1*x(:,1)/r1 + a2*(x(:,1)-xn(1))/r2 + a3*xn(1)/rn
        // f(:,5) = a1*x(:,2)/r1 + a2*(x(:,2)-xn(2))/r2 + a3*xn(2)/rn
        // f(:,6) = a1*x(:,3)/r1 + a2*(x(:,3)-xn(3))/r2 + a3*xn(3)/rn

        for(i=0; i<n; i++)
        {
            tmp1 = this.a1[i]/this.r1[i];
            tmp2 = this.a2[i]/this.r2[i];
            tmp3 = this.a3[i]/this.rn[i];

            f[i][0] = x[i][3];
            f[i][1] = x[i][4];
            f[i][2] = x[i][5];

            f[i][3] = tmp1*x[i][0] + tmp2*(x[i][0]-this.xn[0]) + tmp3*this.xn[0];
            f[i][4] = tmp1*x[i][1] + tmp2*(x[i][1]-this.xn[1]) + tmp3*this.xn[1];
            f[i][5] = tmp1*x[i][2] + tmp2*(x[i][2]-this.xn[2]) + tmp3*this.xn[2];
        }
    }

    /**
     * Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
     */
    this.diffq1 = function(x,f)
    {
        var r21;
        var r1;
        var a1;
        var a2;

        r21 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        r1  = Math.sqrt(r21);
        a1 = -this.m1 / (r21 + this.eps1);
        a2 = -this.m2 / (r21 + this.eps2);
        a1 = a1 + a2;

        f[0] = x[3];
        f[1] = x[4];
        f[2] = x[5];
        f[3] = a1 * x[0] / r1;
        f[4] = a1 * x[1] / r1;
        f[5] = a1 * x[2] / r1;
    }
}
// end SPMModel

// begin NBI model

    var isInit = false;
    var nnn = 10000;
    var nnm1 = nnn-1;
    var rad;
    var rho_halo, mass_halo;
    var rho_disk, mass_disk;
    var rho_bulge, mass_bulge;
    var rho_total, mass_total;
    var masses, radius, density;
    var vr2, vr, new_vr2, new_vr;
    var acceleration, acceleration_particle;
    var new_mass, new_rho, phi;

    var rrout1 = 0;
    var rrout2 = 0;

    var rs_internal = 10.0;
    var rs_internal2 = rs_internal*rs_internal;
    var rs_internal3 = rs_internal*rs_internal*rs_internal;
    var rmax_scale = 100.0;
    var sqrtpi = Math.sqrt(Math.PI);

    var pscale;
    var tmpdfind;

    var lnl;
/**
 * This ForceModel represents an n-body halo/disk/bulge potential that is sampled
 * and then interpolated to calculate the force.
 *
 */
function NBIModel()
{
    // force parameters
    this.params = null;
    this.m1=0;
    this.m2=0;
    this.m3=0;
    this.mm1rs2;
    this.mm2rs2;
    this. mm3rs2;
    this.eps1;
    this.eps2;


    this.pn = 0;
    this.pn1 = 0;
    this.pn2 = 0;


    // temporary storage arrays dependent on n
    this.r22;
    this.r21;
    this.r1;
    this.r2;
    this.rn;
    this.a1;
    this.a2;
    this.a3;
    this.ival11 = null;
    this.ival22 = null;
    this.ivaln = null;
    this.df_force11 = null;
    this.df_force22 = null;
    this.df_forcen = null;
    this.c3n = null;

    // temporary storage independent of n
    this.xn = new Array(6);
    this.xp = 0.0;
    this.yp = 0.0;
    this.zp = 0.0;
    this.vxp = 0.0;
    this.vyp = 0.0;
    this.vzp = 0.0;


    /**
     * Sets the simulation parameters to this model.
     */
    this.setParameters = function(p)
    {
        this.params = p;

        if(this.params != null)
        {
            this.m1 = this.params.mass1;
            this.m2 = this.params.mass2;
            this.m3 = this.params.mass2;
            this.mm1rs2 = -this.m1*rs_internal2;
            this.mm2rs2 = -this.m2*rs_internal2;
            this.mm3rs2 = -this.m3*rs_internal2;
            this.eps1 = this.params.epsilon1 * this.params.epsilon1;
            this.eps2 = this.params.epsilon2 * this.params.epsilon2;
            this.pn = this.params.n;
            this.pn1 = this.params.n1;
            this.pn2 = this.params.n2;
            this.rrout1 = this.params.rout1;
            this.rrout2 = this.params.rout2;
            //if(rad != null)
            //{
                //initDistribution();
            //}
        }

    }


    /** 
     * Initialize temporary storage.
     */
    this.initVars = function(n)
    {
        this.r22 = new Array(n);
        this.r21 = new Array(n);
        this.r1 = new Array(n);
        this.r2 = new Array(n);
        this.rn = new Array(n);
        this.a1 = new Array(n);
        this.a2 = new Array(n);
        this.a3 = new Array(n);
 
        this.ival11 = new Array(n);
        this.ival22 = new Array(n);
        this.ivaln = new Array(n);

        this.df_force11 = new Array(n);
        this.df_force22 = new Array(n);
        this.df_forcen = new Array(n);

        this.c3n = new Array(n);
        for(var i=0; i<n; i++)
        {
        this.r22[i]=0;
        this.r21[i]=0;
        this.r1[i]=0;
        this.r2[i]=0;
        this.rn[i]=0;
        this.a1[i]=0;
        this.a2[i]=0;
        this.a3[i]=0;
 
        this.ival11[i]=0;
        this.ival22[i]=0;
        this.ivaln[i]=0;

        this.df_force11[i]=0;
        this.df_force22[i]=0;
        this.df_forcen[i]=0;

        this.c3n[i]=0;
        }
        this.initDistribution();
    }

    /** 
     * Cleanup temporary storage.
     */
    this.deallocateVars = function()
    {
        this.r22 = null;
        this.r21 = null;
        this.r1 = null;
        this.r2 = null;
        this.rn = null;
        this.a1 = null;
        this.a2 = null;
        this.a3 = null;
 
        this.ival11 = null;
        this.ival22 = null;
        this.ivaln = null;

        this.df_force11 = null;
        this.df_force22 = null;
        this.df_forcen = null;

        this.c3n = null;
    }

    /**
     * Build an n-body halo/disk/bulge
     */
    this.initDistribution = function()
    {
        if(this.isInit) return;
        this.isInit = true;

        var rmax;
        var mold, dmold, mtotal, mtot;
        var rscale;
        var dx, x;
        var alphahalo, qhalo, gammahalo, mhalo, rchalo, rhalo, epsilon_halo;
        var zdisk, zdiskmass, hdisk, zdiskmax;
        var hbulge, mbulge;
        var rho_tmp;
        var G, factor;
        var r, m,  pi;
        var p1, rd, rho_local;
        var p, rr, dr, rh, dp, mnew, dm;
        var acc_merge, rad_merge, acc;
        var xmax;

        var j, nmax, k, nmerge, ntotal, jj;

        rad = new Array(nnn);
        rho_halo = new Array(nnn);
        mass_halo = new Array(nnn);
        rho_disk = new Array(nnn);
        mass_disk = new Array(nnn);
        rho_bulge = new Array(nnn);
        mass_bulge = new Array(nnn);
        rho_total = new Array(nnn);
        mass_total = new Array(nnn);
        masses = new Array(nnn);
        radius = new Array(nnn);
        density = new Array(nnn);
        vr2 = new Array(nnn);
        vr = new Array(nnn);
        new_vr2 = new Array(nnn);
        new_vr = new Array(nnn);
        acceleration = new Array(nnn);
        acceleration_particle = new Array(nnn);
        new_mass = new Array(nnn);
        new_rho = new Array(nnn);
        phi = new Array(nnn);

        // set the constant for dynamical friction
        lnl = 0.00;

        // set up the parameters for the halo
        mhalo = 5.8;
        rhalo = 10.0;
        rchalo = 10.0;
        gammahalo = 1.0;
        epsilon_halo = 0.4;
        pi = Math.PI;

        // derive additional constants
        qhalo = gammahalo / rchalo;
        alphahalo = 1.0 / ( 1.0 - sqrtpi * qhalo * Math.exp(qhalo*qhalo) * (1.0 - this.erf(qhalo)) );

        // set the integration limits and zero integration constants
        rmax = 20;
        nmax = 2000;
        dr = rmax / nmax;
        mold = 0;

        rscale = 5;
        ntotal = nnn;

        // set the limits for integration, and zero integration constants
        k = nmax / 2;
        dx = 1.0  / k;
        x = 0.0;
        dmold = 0.0;
        mtot = 0.0;
        m = 0.0;
        G = 1;

        /////
        // set the fundamental disk parameters
        zdisk = 0.2;
        zdiskmax = 3.5;
        hdisk = 1.0;


        /////
        // set the fundamental bulge parameters
        hbulge = 0.2;
        mbulge = 0.3333;


        ///////////////////////////////////////////////3
        ///// set up the radius array
        for(j = 0; j<nmax; j++)
        {
          x = x + dx;
          rad[j]= x * rchalo;
        }

        ///////////////////////////////////////////////3
        /////

        dr = rad[1] - rad[0];
        dx = dr / rchalo;

        for(j =0; j<nmax; j++)
        {
            // set the position
            r = rad[j];
            x = r / rchalo;

            // calculate the local rho based
            rho_tmp = alphahalo / (2*sqrtpi*sqrtpi*sqrtpi ) * (Math.exp(-x*x) / (x*x + qhalo*qhalo));

            // renormalize for the new halo size
            rho_tmp = rho_tmp / ( rchalo * rchalo * rchalo);

            // calculate mass in local shell, and update total mass
            dm = rho_tmp * 4 * pi * r * r *dr;
            mtot = mtot + dm;

            // store values in an array
            rho_halo[j] = rho_tmp * mhalo;
            mass_halo[j] = mtot * mhalo;
        }

        /////
        // now calculate the potential
        for(j = 0; j<nmax; j++)
        {
            r = rad[j];
            m = mass_halo[j];
            p1 = -G * m / r;
            phi[j] = p1;
        }


        ///////////////////////////////////////////////3
        // disk model

        /////
        // loop over the distribution

        for(j=0; j<nmax; j++)
        {
            // set the radius
            rd = rad[j];

            // find the local density in the disk
            rho_local  = Math.exp(-rd/hdisk)/ (8*pi*hdisk*hdisk);
            rho_disk[j] = rho_local;

            // find the mass in the spherical shell
            mnew = 4 * pi * rho_local * rd *rd * dr;

            mass_disk[j] = mnew + mold;
            mold = mass_disk[j];
        }

        ///////////////////////////////////////////////3
        // bulge model

        /////
        // loop over the distribution
        mold = 0.0;
        for(j=0; j<nmax; j++)
        {
            // set the radius
            rd = rad[j];

            // find the local density in the disk
            rho_local  = Math.exp(-rd*rd/(hbulge*hbulge));
            rho_bulge[j] = rho_local;

            // find the mass in the spherical shell
            mnew = 4 * pi * rho_local * rd *rd * dr;

            mass_bulge[j] = mnew + mold;
            mold = mass_bulge[j];
        }

        // renormalize distribution
        factor = mbulge / mass_bulge[nmax-1];
        for(j=0; j<nmax; j++)
        {
            mass_bulge[j] = mass_bulge[j] * factor;
            rho_bulge[j]  = rho_bulge[j]  * factor;
        }


        dr = rad[1] - rad[0];

        //////////////////////////////////////////////////
        j = 0;
        mass_total[j]=  (mass_halo[j] + mass_disk[j] + mass_bulge[j]);
        r = rad[j];
        rho_total[j] = mass_total[j] /  (4.0/3.0 * pi * r * r * r);
        dr = rad[1] - rad[0];

        for(j = 1; j<nmax; j++)
        {
            r = rad[j];
            mass_total[j]=  (mass_halo[j] + mass_disk[j] + mass_bulge[j]);

            dm = mass_total[j] - mass_total[j-1];
            rho_total[j] = dm / (4 * pi * r * r * dr);
        }

        //////////////////////////////////////////////
        // find the velocity dispersion v_r**2

        masses = mass_total;
        radius = rad;
        density = rho_total;


        for(j=0; j<nmax; j++)
        {
            p = 0.0;
            rr = radius[j];
            dr = radius[nmax-1] / nmax;

            for(jj = j; jj<nmax; jj++)
            {
                m  = masses[jj];
                rh = density[jj];
                rr = rr + dr;

                dp = rh * G * m / (rr*rr) * dr;
                p = p + dp;
            }

            vr2[j] = 1/density[j] * p;
            vr[j] = Math.sqrt(vr2[j]);
        }

        //////////////////////////////////////////////
        // find the velocity dispersion v_r**2
        masses = mass_total;
        radius = rad;
        density = rho_total;

        for(j=0; j<nmax; j++)
        {
            p = 0.0;
            rr = radius[j];
            dr = radius[nmax-1] / nmax;
            for(jj = j; jj<nmax; jj++)
            {
                m  = masses[jj];
                rh = density[jj];
                rr = rr + dr;
                dp = rh * G * m / (rr*rr) * dr;
                p = p + dp;
            }

            vr2[j] = 1/density[j] * p;
            vr[j] = Math.sqrt(vr2[j]);
        }

        //////////////////////////////////////////////
        // find the accelerations felt by the particles and center of mass
        masses = mass_total;
        radius = rad;
        density = rho_total;

        for(j=0; j<nmax; j++)
        {
            rr = radius[j];
            m  = masses[j];
            acceleration[j] = G * m / (rr*rr);
            acceleration_particle[j] = acceleration[j];
        }

        //acceleration_particle = acceleration;
        nmerge = 50;
        acc_merge = acceleration[nmerge];
        rad_merge = rad[nmerge];

        for(j=0; j<nmerge; j++)
        {
            rr = radius[j];
            m  = masses[j];

            // smoothed acceleration
            acc = G * m / (rr*rr + .1* (rad_merge -rr));
            acceleration_particle[j] = acc;
        }

        //////////////////////////////////////////////
        // rederive the masses from the new particle acceleration
        radius = rad;
        dr = rad[1] - rad[0];

        // find the accelerations felt by the particles and center of mass
        radius = rad;

        for(j = 0; j<nmax; j++)
        {
            rr = radius[j];
            new_mass[j] = rr*rr * acceleration_particle[j] / G;
            new_rho[j]  = new_mass[j] / (4 * pi * rr * rr * dr);
        }


        //////////////////////////////////////////////
        // find the velocity dispersion v_r**2 using the new density and masses

        masses = new_mass;
        radius = rad;
        density = new_rho;


        for(j=0; j<nmax; j++)
        {
            p = 0.0;
            rr = radius[j];
            dr = radius[nmax-1] / nmax;

            for(jj=j; jj<nmax; jj++)
            {
                m  = masses[jj];
                rh = density[jj];
                rr = rr + dr;
                dp = rh * G * m / (rr*rr) * dr;
                p = p + dp;
            }

            new_vr2[j] = 1/density[j] * p;
            new_vr[j] = Math.sqrt(new_vr2[j]);
        }


        //////////////////////////////////////////////
        // extend the values to large rmax
        for(j=nmax; j<ntotal; j++)
        {
            mass_total[j] = mass_total[nmax-1];
            mass_halo[j] = mass_halo[nmax-1];
            mass_disk[j] = mass_disk[nmax-1];
            mass_bulge[j] = mass_bulge[nmax-1];
            new_mass[j] = new_mass[nmax-1];

            rho_total[j] = 0.0;
            new_rho[j]   = 0.0;

            vr[j]      = 1e-6;
            vr2[j]     = 1e-6;
            new_vr[j]  = 1e-6;
            new_vr2[j] = 1e-6;

            m = mass_total[nmax-1];
            rr = rad[nmax-1] + dr*(j - nmax);
            rad[j] = rr;
            acc = G * m / (rr*rr);
            acceleration_particle[j] = acc;
            acceleration[j] = acc;
        }

        //////////////////////////////////////////////
        // normalize to the unit mass

        for(j=0; j<ntotal; j++)
        {
            mass_total[j]  = mass_total[j] / 7.13333;
            mass_halo[j]   = mass_halo[j]  / 7.13333;
            mass_disk[j]   = mass_disk[j]  / 7.13333;
            mass_bulge[j]  = mass_bulge[j] / 7.13333;
            new_mass[j]    = new_mass[j]   / 7.13333;

            rho_total[j]   = rho_total[j]  / 7.13333;
            new_rho[j]     = new_rho[j]    / 7.13333;

            vr[j]          = vr[j]      / 7.13333;
            vr2[j]         = vr2[j]     / 7.13333;
            new_vr[j]      = new_vr[j]  / 7.13333;
            new_vr2[j]     = new_vr2[j] / 7.13333;

            rad[j]         = rad[j];

            acceleration_particle[j] = acceleration_particle[j] / 7.13333;
            acceleration[j]          = acceleration[j]  / 7.13333;
        }

        pscale = 1.0;

        tmpdfind = pscale*rs_internal/rmax_scale*nnn;
    } // end initDistribution

    /**
     * Determine the index to use for interpolating the force.
     */
    this.dfIndex = function( rin, rs)
    {
        var ival = 0;
        //double tmp = rin*pscale*rs_internal/rmax_scale;
        var tmp = rin*tmpdfind;
        //tmp = tmp*nnn + 1;
        ival = Math.floor(tmp);
        ival = Math.min(ival, nnn-1);
        //ival = min(int(  (rin * pscale * rs_internal/ rmax_scale) * nnn + 1), nnn)
/*
        ival--; // java is 0-based
        if(ival < 0)
        {
            ival = 0;
        }
*/

        return ival;
    }
   
    /**
     * Computes the error function.  Based on code from Numerical Recipies and the following URL: 
     * http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html
     * fractional error in math formula less than 1.2 * 10 ^ -7.
     * although subject to catastrophic cancellation when z in very close to 0
     * from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
     */
    this.erf = function(z) 
    {
        var t = 1.0 / (1.0 + 0.5 * Math.abs(z));

        // use Horner's method
        var ans = 1 - t * Math.exp( -z*z   -   1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 + 
                                            t * ( 0.09678418 + 
                                            t * (-0.18628806 + 
                                            t * ( 0.27886807 + 
                                            t * (-1.13520398 + 
                                            t * ( 1.48851587 + 
                                            t * (-0.82215223 + 
                                            t * ( 0.17087277))))))))));
        if (z >= 0) return  ans;
        else        return -ans;
    }

    /**
     * Compute the circular velocity for a particle at a distance r from the specified mass.
     * The rout scale of the disk and softening length, eps, are provided.
     */
    this.circularVelocity = function(mass, r, rout, eps)
    {
        var ftotal = 0;

        var ival = this.dfIndex(r, rout);

        ftotal = mass * acceleration_particle[ival] * rs_internal * rs_internal;

        var v = Math.sqrt(ftotal*r);

        return v;
    }

    /**
     * For the given particle positions and velocities, calculate the accelerations.
     */
    this.diffeq = function(x,f)
    {
//Prof.startCall("diffeq");
        var tmp1, tmp2, tmp3, tmp4, tmp5;
        tmp1=tmp2=tmp3=0;
        var n = x.length;
        var i = 0;

        var df_sigma = 0;
        var df_rho = 0;
        var c1 = 0;
        var c2 = 0;
        var xvalue = 0;
        var v1 = 0;
        var v21 = 0;

        for(i=0; i<6; i++)
        {
            this.xn[i] = x[n-1][i];
        }

        xp  = this.xn[0]; yp  = this.xn[1]; zp  = this.xn[2];
        vxp = this.xn[3]; vyp = this.xn[4]; vzp = this.xn[5];

        tmp4 = xp*xp+yp*yp+zp*zp;

        // distance between the two galaxies - the tidal force
        tmp5 = Math.sqrt(tmp4);

var iv1 = 0;
var iv2 = 0;
var ivn = 0;
var df1 = 0;
var df2 = 0;
var dfn = 0;
var tx = null;
var r22d = 0;
var r21d = 0;
var aa1 = 0;
var aa2 = 0;
var aa3 = 0;




        // get the forces, sigma and rho, and rescale them
        //ivaln[0] = dfIndex(rn[0],rrout2);
        ivn = this.dfIndex(tmp5,rrout2);

        df_sigma  =  new_vr2[ivn] * rs_internal2;
        df_rho    =  new_rho[ivn] * rs_internal3;

        // df
        v21 = vxp*vxp + vyp*vyp + vzp*vzp;
        v1  = Math.sqrt(v21);

        xvalue = v1 / df_sigma;
        c1 = this.erf(xvalue) - 2.0 * xvalue / sqrtpi * Math.exp(-xvalue*xvalue);

        // df formula with G=1
        c2  = 4.0 * Math.PI * this.m2 * lnl / v21;
        for(i=0; i<n; i++)
        {
            this.c3n[i] = 0;
        }


        for(i=this.pn1; i<n; i++)
        {
            this.c3n[i]  = c1 * c2 * df_rho;
        }

var tf = null;
var tv1 = 1.0/v1;
var c3tv = 0;
var dv = 0;
        ivn  = this.dfIndex(tmp5, rrout2);

        for(i=0; i<n; i++)
        {
            // distance between the companion and the particle
            tx = x[i];

            //tmp1 = (x[i][0]-xp);
            //tmp2 = (x[i][1]-yp);
            //tmp3 = (x[i][2]-zp);
            tmp1 = (tx[0]-xp);
            tmp2 = (tx[1]-yp);
            tmp3 = (tx[2]-zp);
            //r22[i] = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
            r22d = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

            // distance between the main galaxy and the particle
            //tmp1 = x[i][0];
            //tmp2 = x[i][1];
            //tmp3 = x[i][2];
            tmp1 = tx[0];
            tmp2 = tx[1];
            tmp3 = tx[2];
            //r21[i] = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
            r21d = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;

            //r2[i] = Math.sqrt(r22[i]);
            //r1[i] =  Math.sqrt(r21[i]);
            //r2[i] = Math.sqrt(r22d);
            //r1[i] =  Math.sqrt(r21d);
            //rn[i] = tmp5;

            r22d = Math.sqrt(r22d);
            r21d = Math.sqrt(r21d);
            //r2[i] = r22d;
            //r1[i] = r21d;
            //rn[i] = tmp5;


            //ival11[i]  = dfIndex(r1[i], rrout1);
            //ival22[i] = dfIndex(r2[i], rrout2);
            //ivaln[i]  = dfIndex(rn[i], rrout2);
            //iv1  = dfIndex(r1[i], rrout1);
            //iv2 = dfIndex(r2[i], rrout2);
            //ivn  = dfIndex(rn[i], rrout2);
            //iv1  = dfIndex(r21d, rrout1);
            //iv2 = dfIndex(r22d, rrout2);
            //ivn  = dfIndex(tmp5, rrout2);

            dv = r21d*tmpdfind;
            iv1 = Math.round(dv); if(iv1>nnm1)iv1=nnm1;
            dv = r22d*tmpdfind;
            iv2 = Math.round(dv); if(iv2>nnm1)iv2=nnm1;

            //df_force11[i] = acceleration_particle[ival11[i]] * rs_internal2;
            //df_force22[i] = acceleration_particle[ival22[i]] * rs_internal2;
            //df_forcen[i]  = acceleration_particle[ivaln[i]]  * rs_internal2;
            //df1 = acceleration_particle[iv1] * rs_internal2;
            //df2 = acceleration_particle[iv2] * rs_internal2;
            //dfn  = acceleration_particle[ivn]  * rs_internal2;

            //a1[i] = -m1 * df_force11[i];
            //a2[i] = -m2 * df_force22[i];
            //a3[i] = -m3 * df_forcen[i];
            //a1[i] = -m1 * df1;
            //a2[i] = -m2 * df2;
            //a3[i] = -m3 * dfn;
            //a1[i] = acceleration_particle[iv1]*mm1rs2;
            //a2[i] = acceleration_particle[iv2]*mm2rs2;
            //a3[i] = acceleration_particle[ivn]*mm3rs2;
            aa1 = acceleration_particle[iv1]*this.mm1rs2;
            aa2 = acceleration_particle[iv2]*this.mm2rs2;
            aa3 = acceleration_particle[ivn]*this.mm3rs2;

            //tmp1 = a1[i]/r1[i];
            //tmp2 = a2[i]/r2[i];
            //tmp3 = a3[i]/rn[i];
            tmp1 = aa1/r21d;
            tmp2 = aa2/r22d;
            tmp3 = aa3/tmp5;
            c3tv = this.c3n[i]*tv1;

            //f[i][0] = x[i][3];
            //f[i][1] = x[i][4];
            //f[i][2] = x[i][5];

            //f[i][3] = tmp1*x[i][0] + tmp2*(x[i][0]-xn[0]) + tmp3*xn[0] - c3n[i] * xn[3] / v1;
            //f[i][4] = tmp1*x[i][1] + tmp2*(x[i][1]-xn[1]) + tmp3*xn[1] - c3n[i] * xn[4] / v1;
            //f[i][5] = tmp1*x[i][2] + tmp2*(x[i][2]-xn[2]) + tmp3*xn[2] - c3n[i] * xn[5] / v1;

            tf = f[i];
            //tx = x[i];
            tf[0] = tx[3];
            tf[1] = tx[4];
            tf[2] = tx[5];

            tf[3] = tmp1*tx[0] + tmp2*(tx[0]-this.xn[0]) + tmp3*this.xn[0] - this.xn[3] * c3tv;
            tf[4] = tmp1*tx[1] + tmp2*(tx[1]-this.xn[1]) + tmp3*this.xn[1] - this.xn[4] * c3tv;
            tf[5] = tmp1*tx[2] + tmp2*(tx[2]-this.xn[2]) + tmp3*this.xn[2] - this.xn[5] * c3tv;
        }

// redo acceleration for last particle
r22d=1.0;
tmp2 = aa2/r22d;
            tf[3] = tmp1*tx[0] + tmp2*(tx[0]-this.xn[0]) + tmp3*this.xn[0] - this.xn[3] * c3tv;
            tf[4] = tmp1*tx[1] + tmp2*(tx[1]-this.xn[1]) + tmp3*this.xn[1] - this.xn[4] * c3tv;
            tf[5] = tmp1*tx[2] + tmp2*(tx[2]-this.xn[2]) + tmp3*this.xn[2] - this.xn[5] * c3tv;


//Prof.endCall("diffeq");
    }

    /**
     * Calculate the acceleration felt by the secondary galaxy in orbit around the primary.
     */
    this.diffq1 = function(x,f)
    {
        var r21, r1, a1, a2, at;
        var c1, c2, c3, v21, v1, xvalue;

        var ival, ival2;
        var df_force1, df_force2;
        var df_sigma, df_rho;

        var i;
        var rr, rlocal, ivalue, dr, mmm, dm;

        r21 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        r1  = Math.sqrt(r21);

        // get the index for the calculations
        ival  =  this.dfIndex(r1, rrout1   );
        ival2 =  this.dfIndex(r1, rrout2   );

        // get the forces, sigma and rho, and rescale them
        df_force1 = acceleration_particle[ival] * rs_internal * rs_internal;
        df_force2 = acceleration_particle[ival2]* rs_internal * rs_internal;

        df_sigma  = new_vr2[ival] * rs_internal * rs_internal;
        df_rho    = new_rho[ival] * ( rs_internal * rs_internal * rs_internal );

        // interpolated forces
        a1 = -this.m1 * df_force1;
        a2 = -this.m2 * df_force2;
        at = a1 + a2;

        // df
        v21 = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
        v1  = Math.sqrt(v21);

        xvalue = v1 / df_sigma;
        c1 = this.erf(xvalue) - 2.0 * xvalue / sqrtpi * Math.exp(-xvalue*xvalue);

        // df formula with G=1
        c2 = -4.0 * Math.PI * this.m2 * lnl / v21;
        c3 = c1 * c2 * df_rho;

        f[0] = x[3];
        f[1] = x[4];
        f[2] = x[5];

        f[3] = at * x[0]/r1 - c3 * x[3]/v1;
        f[4] = at * x[1]/r1 - c3 * x[4]/v1;
        f[5] = at * x[2]/r1 - c3 * x[5]/v1;
    }
}
// end NBI model

// begin integrator
function Integrator()
{
    var od6 = 1.0/6.0;
    var od3 = 1.0/3.0;

    this.force;

    /** temporary array used in integration */
    this.x;

    /** temporary array used in integration */
    this.xe;

    /** internal array to hold force,acceleration calculations */
    this.f;

    this.setForce = function(forceIn)
    {
        this.force = forceIn;
    } 

    /**
     * Initialize temporary storage arrays
     */
    this.initRKVar = function(n)
    {
        this.x  = new Array(n+1);
        this.xe  = new Array(n+1);
        this.f  = new Array(n+1);
      
        for(var i=0; i<=n; i++)
        {
            this.x[i] = new Array(6); 
            this.xe[i] = new Array(6); 
            this.f[i] = new Array(6); 
for(var j=0; j<6; j++)
{
this.x[i][j]=0;
this.xe[i][j]=0;
this.f[i][j]=0;
}
        }

        if(this.force != null)
        {
            this.force.initVars(n+1);
        }
    }

    /**
     * Cleanup temporary storage.
     */
    this.deallocateRKVar = function()
    {
        this.x = null;
        this.xe = null;
        this.f = null;

        if(this.force != null)
        {
            this.force.deallocateVars();
        }
    }

    /**
     * Using the provided ForceModel, compute the updated positions and velocities
     * by taking a step of size h.  These are computed for all particles.
     */
    this.rk4 = function(x0, xout, h)
    {
        var n = x0.length;
        var i=0;
        // pre-compute products
        var hod6 = h*od6;
        var h0p5 = h*0.5;
        var hod3 = h*od3;
        // we unroll the second dimension of 6 elements to speed computation       
//Prof.startCall("rk4_1");
        for(i=0;i<n;i++)
        {
            this.x[i][0] = x0[i][0];
            this.x[i][1] = x0[i][1];
            this.x[i][2] = x0[i][2];
            this.x[i][3] = x0[i][3];
            this.x[i][4] = x0[i][4];
            this.x[i][5] = x0[i][5];
        }
//Prof.endCall("rk4_1");
          
        this.force.diffeq(this.x,this.f);

//Prof.startCall("rk4_2");
        for(i=0;i<n;i++)
        {
            this.xe[i][0] = x0[i][0] + this.f[i][0] * hod6;
            this.x[i][0]  = x0[i][0] + this.f[i][0] * h0p5;
            this.xe[i][1] = x0[i][1] + this.f[i][1] * hod6;
            this.x[i][1]  = x0[i][1] + this.f[i][1] * h0p5;
            this.xe[i][2] = x0[i][2] + this.f[i][2] * hod6;
            this.x[i][2]  = x0[i][2] + this.f[i][2] * h0p5;
            this.xe[i][3] = x0[i][3] + this.f[i][3] * hod6;
            this.x[i][3]  = x0[i][3] + this.f[i][3] * h0p5;
            this.xe[i][4] = x0[i][4] + this.f[i][4] * hod6;
            this.x[i][4]  = x0[i][4] + this.f[i][4] * h0p5;
            this.xe[i][5] = x0[i][5] + this.f[i][5] * hod6;
            this.x[i][5]  = x0[i][5] + this.f[i][5] * h0p5;
        }
//Prof.endCall("rk4_2");

        this.force.diffeq(this.x,this.f);

//Prof.startCall("rk4_3");
        for(i=0;i<n;i++)
        {
            this.xe[i][0] = this.xe[i][0] + this.f[i][0] * hod3;
            this.x[i][0]  = x0[i][0] + this.f[i][0] * h0p5;
            this.xe[i][1] = this.xe[i][1] + this.f[i][1] * hod3;
            this.x[i][1]  = x0[i][1] + this.f[i][1] * h0p5;
            this.xe[i][2] = this.xe[i][2] + this.f[i][2] * hod3;
            this.x[i][2]  = x0[i][2] + this.f[i][2] * h0p5;
            this.xe[i][3] = this.xe[i][3] + this.f[i][3] * hod3;
            this.x[i][3]  = x0[i][3] + this.f[i][3] * h0p5;
            this.xe[i][4] = this.xe[i][4] + this.f[i][4] * hod3;
            this.x[i][4]  = x0[i][4] + this.f[i][4] * h0p5;
            this.xe[i][5] = this.xe[i][5] + this.f[i][5] * hod3;
            this.x[i][5]  = x0[i][5] + this.f[i][5] * h0p5;
        }
//Prof.endCall("rk4_3");

        this.force.diffeq(this.x,this.f);

//Prof.startCall("rk4_4");
        for(i=0;i<n;i++)
        {
            this.xe[i][0] = this.xe[i][0] + this.f[i][0] * hod3;
            this.x[i][0]  = x0[i][0] + h * this.f[i][0];
            this.xe[i][1] = this.xe[i][1] + this.f[i][1] * hod3;
            this.x[i][1]  = x0[i][1] + h * this.f[i][1];
            this.xe[i][2] = this.xe[i][2] + this.f[i][2] * hod3;
            this.x[i][2]  = x0[i][2] + h * this.f[i][2];
            this.xe[i][3] = this.xe[i][3] + this.f[i][3] * hod3;
            this.x[i][3]  = x0[i][3] + h * this.f[i][3];
            this.xe[i][4] = this.xe[i][4] + this.f[i][4] * hod3;
            this.x[i][4]  = x0[i][4] + h * this.f[i][4];
            this.xe[i][5] = this.xe[i][5] + this.f[i][5] * hod3;
            this.x[i][5]  = x0[i][5] + h * this.f[i][5];
        }
//Prof.endCall("rk4_4");

        this.force.diffeq(this.x,this.f);

//Prof.startCall("rk4_5");
        // combine last sub-step with copying to output array
        for(i=0;i<n;i++)
        {
            this.xe[i][0] = this.xe[i][0] + this.f[i][0] * hod6;
            this.xe[i][1] = this.xe[i][1] + this.f[i][1] * hod6;
            this.xe[i][2] = this.xe[i][2] + this.f[i][2] * hod6;
            this.xe[i][3] = this.xe[i][3] + this.f[i][3] * hod6;
            this.xe[i][4] = this.xe[i][4] + this.f[i][4] * hod6;
            this.xe[i][5] = this.xe[i][5] + this.f[i][5] * hod6;
              
            xout[i][0] = this.xe[i][0];
            xout[i][1] = this.xe[i][1];
            xout[i][2] = this.xe[i][2];
            xout[i][3] = this.xe[i][3];
            xout[i][4] = this.xe[i][4];
            xout[i][5] = this.xe[i][5];
        }
//Prof.endCall("rk4_5");

        return;
    } // end rk4

    /**
     * Using the provided ForceModel, compute the updated position and velocity of
     * the secondary galaxy by taking a step of size h.
     */
    this.rk41 = function(xx0, xxe, h)
    {
        var x = new Array(7);
        var f = new Array(7);
        var i=0;
        var n=6;

        // pre-compute products
        var hod6 = h*od6;
        var h0p5 = h*0.5;
        var hod3 = h*od3;

        //for(i=0;i<n;i++)
        //{
                //x[i] = xx0[i];
        //}

        x[0] = xx0[0];
        x[1] = xx0[1];
        x[2] = xx0[2];
        x[3] = xx0[3];
        x[4] = xx0[4];
        x[5] = xx0[5];
        x[6] = 0;
        f[0] = f[1] = f[2] = f[3] = f[4] = f[5] = f[6] = 0;
        this.force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
            //xxe[i] = xx0[i] + h * f[i] * od6;
            //x[i]   = xx0[i] + h * f[i] * 0.5;
            xxe[0] = xx0[0] + f[0] * hod6;
            x[0]   = xx0[0] + f[0] * h0p5;
            xxe[1] = xx0[1] + f[1] * hod6;
            x[1]   = xx0[1] + f[1] * h0p5;
            xxe[2] = xx0[2] + f[2] * hod6;
            x[2]   = xx0[2] + f[2] * h0p5;
            xxe[3] = xx0[3] + f[3] * hod6;
            x[3]   = xx0[3] + f[3] * h0p5;
            xxe[4] = xx0[4] + f[4] * hod6;
            x[4]   = xx0[4] + f[4] * h0p5;
            xxe[5] = xx0[5] + f[5] * hod6;
            x[5]   = xx0[5] + f[5] * h0p5;
        //}

        this.force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
                //xxe[i] = xxe[i] + h * f[i] * od3;
                //x[i]   = xx0[i] + h * f[i] *0.5;
                xxe[0] = xxe[0] + f[0] * hod3;
                x[0]   = xx0[0] + f[0] * h0p5;
                xxe[1] = xxe[1] + f[1] * hod3;
                x[1]   = xx0[1] + f[1] * h0p5;
                xxe[2] = xxe[2] + f[2] * hod3;
                x[2]   = xx0[2] + f[2] * h0p5;
                xxe[3] = xxe[3] + f[3] * hod3;
                x[3]   = xx0[3] + f[3] * h0p5;
                xxe[4] = xxe[4] + f[4] * hod3;
                x[4]   = xx0[4] + f[4] * h0p5;
                xxe[5] = xxe[5] + f[5] * hod3;
                x[5]   = xx0[5] + f[5] * h0p5;
        //}

        this.force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
                //xxe[i] = xxe[i] + h * f[i] * od3;
                //x[i]   = xx0[i] + h*f[i];
                xxe[0] = xxe[0] + f[0] * hod3;
                x[0]   = xx0[0] + h*f[0];
                xxe[1] = xxe[1] + f[1] * hod3;
                x[1]   = xx0[1] + h*f[1];
                xxe[2] = xxe[2] + f[2] * hod3;
                x[2]   = xx0[2] + h*f[2];
                xxe[3] = xxe[3] + f[3] * hod3;
                x[3]   = xx0[3] + h*f[3];
                xxe[4] = xxe[4] + f[4] * hod3;
                x[4]   = xx0[4] + h*f[4];
                xxe[5] = xxe[5] + f[5] * hod3;
                x[5]   = xx0[5] + h*f[5];
        //}

        this.force.diffq1(x, f);

        //for(i=0;i<n;i++)
        //{
                //xxe[i] = xxe[i] + h * f[i] * od6;
                xxe[0] = xxe[0] + f[0] * hod6;
                xxe[1] = xxe[1] + f[1] * hod6;
                xxe[2] = xxe[2] + f[2] * hod6;
                xxe[3] = xxe[3] + f[3] * hod6;
                xxe[4] = xxe[4] + f[4] * hod6;
                xxe[5] = xxe[5] + f[5] * hod6;
        //}

        return;
    } // end rk41
}

// end integrator

function parseDouble(str)
{
    return parseFloat(str);
}

function parseStateInfoString(params,str)
{
        var vals = new Array(22);
        var sa = str.split(",");
        var len = Math.min(sa.length,vals.length);
        for(var i=0; i<len; i++)
        {
            vals[i] = parseDouble(sa[i]);
        }

        params.sec_vec[0] = vals[0];
        params.sec_vec[1] = vals[1];
        params.sec_vec[2] = vals[2];
        params.sec_vec[3] = vals[3];
        params.sec_vec[4] = vals[4];
        params.sec_vec[5] = vals[5];

        params.mass1 = vals[6];
        params.mass2 = vals[7];
        params.rout1 = vals[8];
        params.rout2 = vals[9];
        params.phi1 = vals[10];
        params.phi2 = vals[11];
        params.theta1 = vals[12];
        params.theta2 = vals[13];
        params.epsilon1 = vals[14];
        params.epsilon2 = vals[15];
        params.rscale1[0] = vals[16];
        params.rscale1[1] = vals[17];
        params.rscale1[2] = vals[18];
        params.rscale2[0] = vals[19];
        params.rscale2[1] = vals[20];
        params.rscale2[2] = vals[21];

        params.use_sec_vec = true;
}
// begin setuputil

function SetupUtil()
{
    this.forceModel;
    this.integrator;
    this.params;

    this.createForceModel = function(potential_type, apply)
    {
        var force = null;
        switch(potential_type)
        {
            case 0:
            default:
                force = new SPMModel();
                break;
            case 1:
                force = new NBIModel();
                break;
/*
            case 2:
                force = new MONDModel();
                break;
*/
        }

        if(apply)
        {
            this.forceModel = force;
            if(this.integrator != null) this.integrator.force=force;
        }

        return force;
    }

    this.setHelpers = function(fm, intg, p)
    {
        this.forceModel = fm;
        this.params = p;
        this.integrator = intg;
    }

    /**
     * Determine if the caller provides a parameter file or a parameter string
     */
    this.customCollision = function(args)
    {
        this.params.tIsSet = false;

        if(args.length > 0)
        {
            //if(args.length > 1 && args[0].toLowerCase().equals("-f"))
            //{
                // parse the file
                //IOUtil.readParameterFile(params,args[1]);
            //}
            //else
            //{
                // parse the state string
                parseStateInfoString(this.params,args[0]);
                this.params.potential_type=1;
                this.params.h = 0.001;//Parameters.hbase;
                this.params.tstart = -5;
                this.params.tend = 0;
                this.params.time = -5;

                if(args.length > 1)
                { 
                    var t = parseDouble(args[1]);
                    if(t != 0)
                    {
                        this.params.tstart = t;
                        this.params.time = t; 
                        this.params.tIsSet = true;
                    }
                }
            //}
        }

        this.params.n = this.params.n1 + this.params.n2;
        this.params.eps1 = this.params.epsilon1*this.params.epsilon1;
        this.params.eps2 = this.params.epsilon2*this.params.epsilon2;

    }

    /**
     * Initialize the disks and set the initial positions.
     */
    this.createCollision = function()
    {
        this.profile(this.params.rin1, this.params.rout1, this.params.rscale1, 0, this.params.n1, this.params.mass1,
                this.params.eps1, this.params.theta1, this.params.phi1, this.params.opt1, this.params.heat1,
                this.integrator.x);        

        this.profile(this.params.rin2, this.params.rout2, this.params.rscale2, this.params.n1, this.params.n, this.params.mass2,
                this.params.eps2, this.params.theta2, this.params.phi2, this.params.opt2, this.params.heat2,
                this.integrator.x);        

        // determine if we need to calculate tStart
/*
        if( !this.params.tIsSet ) 
        {
            var sec_vec = this.params.sec_vec;
            double rv4min[] = {sec_vec[0],sec_vec[1],sec_vec[2],-sec_vec[3],-sec_vec[4],-sec_vec[5],0.0};
            double tminVals[] = getTStart(rv4min, params.mass1, params.mass2, params.eps1, params.eps2 , params.h,-30.0 ,10.0*params.rout1, params.rout1, params.rout2);

            double tmpT = tminVals[0];
            if ( tmpT < -12.0 ) 
            {
                tmpT = -5;
            }

            if ( Math.abs(tmpT) < params.h)
            {
                tmpT = -5;
            }
            params.tstart = tmpT;
            params.time = tmpT;
            params.tIsSet  = true;
        }
*/
        var mins = null;
        if(this.params.use_sec_vec)
        {
            mins = this.perturberPositionVec(this.params.sec_vec, this.params.mass1, this.params.mass2, 
                                        this.params.eps1, this.params.eps2,
                                        this.params.h, this.params.n, this.params.n1, this.params.time, this.integrator.x);
        }
        else
        {
            mins = this.perturberPosition(this.params.inclination_degree, this.params.omega_degree, 
                                     this.params.rmin, this.params.velocity_factor, this.params.rout1,
                                     this.params.mass1, this.params.mass2, this.params.eps1, this.params.eps2,
                                     this.params.h, this.params.n, this.params.n1, this.params.time, 
                                     this.integrator.x, this.params.sec_vec);
        }

        this.params.tmin = mins[0];
        this.params.rmin = mins[1];
        this.params.vrmin = mins[2];
        this.params.beta = mins[3];
    }

    this.randm = function()
    {
        return Math.random();
        //return rand.nextDouble();
    }

    this.distrb = function(r1, opt, rscale)
    {
    	  var distrb = 0;
    	  
    	  if (opt == 1)
    	  {
    	    distrb = 1.0/r1;
    	  }
    	  else if (opt == 2)
    	  {
    	    distrb = Math.exp(-r1/rscale[0]);
    	  }
    	  else if (opt == 3)
    	  {
    	    distrb = Math.exp(-r1*r1*rscale[0] - rscale[1]*r1 - rscale[2] );
    	  }
    	  
    	  return distrb;
    }//	end function distrb


    this.profile = function(rin, rout, rscale, nstart, ntot, 
                        mass, eps, theta, phi, opt, 
                        heat, x0)
    {
        var stheta,ctheta,sphi,cphi,pi;
        var x3,y3,z3,xv3,yv3,zv3,x2,y2,z2,xv2,yv2,zv2;
        var x,y,z,xv,yv,zv;

        var i, j, n, nprof;

        var rnorm;
        var rp, r, angle, v, p_ring, cp_ring;
        var st, ct, dr, ran, r1, r2, ptot;
        var n_ring;
        var nring, dnring, is, ie, iring, tct;

        var xp, yp, zp;
        var fx, fy, fz, ftotal;
        var tmp;

        var ival;

        n = x0.length;
if(ntot > n)
{
ntot = n;
alert("ntot exceeds n");
}
        r = new Array(n);
        angle = new Array(n); 
        v = new Array(n);

        pi = Math.PI;

        stheta = Math.sin(theta*pi/180.0);    
        ctheta = Math.cos(theta*pi/180.0);    
        sphi   = Math.sin(phi*pi/180.0);    
        cphi   = Math.cos(phi*pi/180.0);    

        // set up the probability distribution for the disk
        nprof = 1000;
        nring = nprof / 10;

        dnring = nprof/nring;
        rp = new Array(nprof);
        n_ring = new Array(nprof);
        p_ring = new Array(nprof);
        cp_ring = new Array(nprof);

        // set the differential sum of the probability funtion into a vector
        rnorm = 0.0;
        dr = (rout - rin)/(nprof);
        for(i=0; i<nprof; i++)
        {
            r1 = i*dr + rin;
            rp[i] = this.distrb(r1, opt, rscale) * r1 * dr * 2.0 * pi;
            rnorm = rnorm + rp[i];
        }

        // normalize the vector
        for(i=0; i<nprof; i++)
        {
            rp[i] /= rnorm;
        }

        // take the fine bins and put them into the selection bins
        tct = 0;
        for(iring = 0; iring < nring; iring++)
        {
        //do iring =  1, nring
          is = (iring) * dnring + 1;
          ie = (iring+1) * dnring ;

          ptot = 0.0;
          for(i=is;i<ie;i++)
          {
          //do i = is, ie
            ptot = ptot + rp[i];
          } //enddo
          p_ring[iring] = ptot;
        } //enddo

    // formulative cumulative distribution function
    //
        cp_ring[0] = p_ring[0];
        for(iring=1;iring<nring;iring++)
        {
        //do iring = 2, nring
          cp_ring[iring] = cp_ring[iring -1] + p_ring[iring];

        } //enddo


    // find the number of particles in each bin
    //
        n_ring = new Array(nprof);
for(i=0; i<nprof; i++)n_ring[i] = 0;
        for(i=nstart-1;i<ntot;i++)
        {

            // find the radial position bin
            ran = this.randm();
            j = 0;
           // while (j < nprof && ran > rp[j]/rnorm)
                //while (j < nprof && ran > rp[j])
            while(j<nring && ran > cp_ring[j])
            {
              j = j + 1;
            } //enddo
            j--;
            if(j<0)
           {
                j=0;
            }
            n_ring[j] = n_ring[j] + 1;
        } //enddo

        tct = 0;
        i = nstart;
        for(iring=0;iring<nring&&i<n;iring++)
        {
        //do iring =  1, nring
          is = (iring) * dnring + 1;
          ie = (iring+1) * dnring;

          r1 = (is)*dr + rin;
          r2 = (ie)*dr + rin;
          for(j=0; j<n_ring[iring]&&i<n; j++)
          {
          //do j = 1, n_ring(iring)
            ran = this.randm();
            r[i] = r1 + ran * (r2 - r1);
            i = i + 1;
          } //enddo
        } //enddo

    //   set the angular positions and orbital velocities
    //
        for(i=nstart;i<ntot;i++)
        {
        //do i=nstart,ntot
          angle[i] = 2.0 * pi * this.randm();
            v[i] = this.forceModel.circularVelocity(mass,r[i],rout,eps);
        }

        // set position and velocity based on the distribution parameters
        for(i= nstart; i<ntot; i++)
        {
            st  =  Math.sin(angle[i]);
            ct  =  Math.cos(angle[i]);

            x   =  ct*r[i];
            y   =  st*r[i];
            z   =  0.0;

            xv  = -v[i]*st;
            yv  =  v[i]*ct;
            zv  =  0.0;

            x2  =   x * ctheta +  z * stheta;
            y2  =   y;
            z2  =  -x * stheta +  z * ctheta;
            xv2 =  xv * ctheta + zv * stheta;
            yv2 =  yv;
            zv2 = -xv * stheta + zv * ctheta;

            x3  =  x2  * cphi -  y2 * sphi;
            y3  =  x2  * sphi +  y2 * cphi;
            z3  =  z2;
            xv3 =  xv2 * cphi - yv2 * sphi;
            yv3 =  xv2 * sphi + yv2 * cphi;
            zv3 =  zv2;
      

            x0[i][0] = x3;
            x0[i][1] = y3;
            x0[i][2] = z3;
            x0[i][3] = xv3  + this.randm()*heat;
            x0[i][4] = yv3  + this.randm()*heat;
            x0[i][5] = zv3  + this.randm()*heat;

        }
    }

    /**
     *
     */
    this.perturberPosition = function(inclinationDegree, omegaDegree, 
                                      rMin, velocityFactor, rout1,
                                      mass1, mass2, eps1, eps2,
                                      h, n, n1, t0, x0, sec_vec)
    {
        var xx0 = new Array(6);

        var omega = Math.toRadians(omegaDegree);
        var inc = Math.toRadians(inclinationDegree);

	var v = Math.sqrt(2.0)*this.forceModel.circularVelocity(mass1+mass2,rMin,rout1,eps1);

        v = -v * velocityFactor;

        //      setup the transformaton based on the matrix in
        //      fundamentals of astrodynamics p. 82 by
        //      bate, mueller, and white (1971)
        //
        xx0[0] = Math.cos(omega) * rMin;
        xx0[1] = Math.sin(omega) * Math.cos(inc) * rMin;
        xx0[2] = Math.sin(omega) * Math.sin(inc) * rMin;

        xx0[3] = -Math.sin(omega) * v;
        xx0[4] =  Math.cos(omega) * Math.cos(inc) * v;
        xx0[5] =  Math.cos(omega) * Math.sin(inc) * v;

        sec_vec[0] = xx0[0];
        sec_vec[1] = xx0[1];
        sec_vec[2] = xx0[2];
        sec_vec[3] = -xx0[3];
        sec_vec[4] = -xx0[4];
        sec_vec[5] = -xx0[5];

        return this.perturberPositionVec(sec_vec, mass1, mass2, eps1, eps2, h, n, n1, t0, x0);
    }

    this.perturberPositionVec = function(xx0, mass1, mass2, 
                                         eps1, eps2, 
                                         h, n, n1, t0, x0)
    {
        var xxe = new Array(6);
        var tcurrent = 0;
        var i = 0;

        // reverse the velocity for backward integration
        xx0[3] = -xx0[3];
        xx0[4] = -xx0[4];
        xx0[5] = -xx0[5];

        var tmin = tcurrent;
        var tmpr = 0;
        var tmpv = 0;
        // avoid multiple calls to sqrt, save it until the end
        var rmin = xx0[0]*xx0[0]+xx0[1]*xx0[1]+xx0[2]*xx0[2];
        var vrmin = xx0[3]*xx0[3]+xx0[4]*xx0[4]+xx0[5]*xx0[5];
        // now move position back to t0 from t=0.0
        while( t0 < tcurrent )
        {
            this.integrator.rk41(xx0, xxe, h);
            xx0[0] = xxe[0];  xx0[1] = xxe[1]; xx0[2] = xxe[2];
            xx0[3] = xxe[3];  xx0[4] = xxe[4]; xx0[5] = xxe[5];

            tcurrent = tcurrent - h;

            tmpr = xx0[0]*xx0[0]+xx0[1]*xx0[1]+xx0[2]*xx0[2];
            if(tmpr < rmin)
            {
                rmin = tmpr;
                vrmin = xx0[3]*xx0[3]+xx0[4]*xx0[4]+xx0[5]*xx0[5];
                tmin = this.tcurrent;
            }
        }

        // reverse the velocity for forward integration
        xx0[3] = -xx0[3];
        xx0[4] = -xx0[4];
        xx0[5] = -xx0[5];
        // now adjust the test particles from the second disk
        // to the proper velocity and positions
        if(this.params.n > this.params.n1)
        {
            for(i=this.params.n1; i<this.params.n; i++)
            {
                x0[i][0] += xx0[0];
                x0[i][1] += xx0[1];
                x0[i][2] += xx0[2];
                x0[i][3] += xx0[3];
                x0[i][4] += xx0[4];
                x0[i][5] += xx0[5];
            }
        }

        // include the perturbing galaxy
        i = this.params.n;
        x0[i][0] = xx0[0];
        x0[i][1] = xx0[1];
        x0[i][2] = xx0[2];
        x0[i][3] = xx0[3];
        x0[i][4] = xx0[4];
        x0[i][5] = xx0[5];

        vrmin = Math.sqrt(vrmin);
        // beta = (m1+m2)/(rmin*rmin + vrmin)
        var beta = (mass1+mass2)/(rmin*vrmin);
        rmin = Math.sqrt(rmin);
        //return new double[]{tmin,rmin,vrmin,beta};
        var out = new Array(4);
        out[0] = tmin;
        out[1] = rmin;
        out[2] = vrmin;
        out[3] = beta;
        return out;
    }

}

// end setuputil

// begin run
function SPAMRun()
{
    this.params;
    this.forceModel;
    this.integrator;
    this.x0;
    this.xout;

    this.init = function()
    {
        this.params = new Parameters();
        this.params.defaultParameters();
        //this.forceModel = new SPMModel();
        this.forceModel = new NBIModel();
        this.integrator = new Integrator();
        this.integrator.setForce(this.forceModel);
    }

    this.initRun = function(args)
    {
//Prof.startCall("initRun(args[])");
        var su = new SetupUtil();

        su.setHelpers(this.forceModel, this.integrator, this.params);
        su.customCollision(args);
        // update the forceModel based upon passed in args
        this.forceModel = su.createForceModel(this.params.potential_type,true);

        this.forceModel.setParameters(this.params);
        this.integrator.initRKVar(this.params.n);

        this.x0  = new Array(this.params.n+1);
        this.xout = new Array(this.params.n+1);
        var len = this.params.n+1;
        for(var i=0; i<len; i++)
        {
            this.x0[i] = new Array(6);
            this.xout[i] = new Array(6);
for(var j=0; j<6; j++)
{
this.x0[i][j]=0;
this.xout[i][j]=0;
}
        }

        su.createCollision();

        this.copyParticles(this.integrator.x,this.x0);

        this.forceModel.setParameters(this.params);
//Prof.endCall("initRun(args[])");
    }

    /**
     * Copies the particles from x1 to x2.
     */
    this.copyParticles = function(x1, x2)
    {
//Prof.startCall("copyParticles");
        var n = x1.length;
        // assuming 6 members

        for(var i=0; i<n; i++)
        {
            x2[i][0] = x1[i][0];
            x2[i][1] = x1[i][1];
            x2[i][2] = x1[i][2];
            x2[i][3] = x1[i][3];
            x2[i][4] = x1[i][4];
            x2[i][5] = x1[i][5];
        }
//Prof.endCall("copyParticles");
    }

    this.takeAStep = function()
    {
//Prof.startCall("takeAStep");
        this.integrator.rk4(this.x0,this.xout,this.params.h);
        this.copyParticles(this.xout,this.x0);
        this.params.time = this.params.time + this.params.h;
//Prof.endCall("takeAStep");
    }

}

// end run
