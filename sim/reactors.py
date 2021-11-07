import numpy as np
import cantera as ct


def st(temp, pressure, mix, gas, target_spcs, end_time, p_of_t=None,
       run_rop=False):
    """ Runs a shock tube simulation, with non-ideal pressure rise, dP/dt,
        incorporated via an input P(t) profile

        :param temp: reactor inlet temperature (Kelvin)
        :type temp: float
        :param pressure: reactor constant pressure (atm)
        :type pressure: float
        :param mix: initial species concentrations at reactor inlet
        :type mix: dct {spc1: conc1, spc2: ...}
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param target_spcs: desired species concentrations
        :type target_spcs: list [spc1, spc2, ...]
        :param end_time: end time for the simulation (s); simulation will end
            at the first of either this time or the max time in p_of_t
        :type end_time: float
        :param p_of_t: specified pressure (any units) vs time (s) profile;
            can be None if no non-ideal effects are being simulated
        :type p_of_t: Numpy array of shape (2, num_exp_timesteps)
        :return target_concs: target_spcs mole fractions at each time
        :rtype: Numpy array of shape (num_target_spcs, num_timesteps)
        :return times: array of time values (s)
        :rtype: Numpy array of shape (num_timesteps,)

        Note: the timesteps of the p_of_t and the simulation results do not have
            to be (and likely will not be) the same length
    """

    def _wall_velocity(time):
        """ Gets the wall velocity at some time based on v_of_t

            :param time: individual time value (s)
            :type time: float
            :return vel: velocity of the wall (m/s)
            :rtype: float
        """
        dvdt = np.gradient(v_of_t[1, :], v_of_t[0, :])
        dxdt = dvdt / wall.area  # only for clarity (since A = 1 m^2)
        vel = np.interp(time, v_of_t[0, :], -1 * dxdt)

        return vel

    pressure = pressure * 101325  # atm to Pa
    gas.TPX = temp, pressure, mix

    if p_of_t is None:
        # Create an array of ones spanning from t=0 to t=end_time
        v_of_t = np.array([[0, end_time], [1, 1]])
    else:
        # Convert P(t) to V(t) assuming isentropic behavior
        gamma = gas.cp_mass / gas.cv_mass  # specific heat ratio
        v_of_t = np.zeros_like(p_of_t)
        v_of_t[0, :] = p_of_t[0, :]  # time values are the same
        v_of_t[1, :] = (p_of_t[1, :] / p_of_t[1, 0]) ** (-1 / gamma)

    # Set up reactor, wall, and network
    reac = ct.IdealGasReactor(gas)
    reac.volume = 1  # m^3 (this is the default, but put here for clarity)
    env = ct.Reservoir(gas)
    wall = ct.Wall(env, reac)
    wall.area = 1  # m^2 (this is the default, but put here for clarity)
    wall.heat_transfer_coeff = 0  # no heat transfer considerations for ST
    wall.set_velocity(_wall_velocity)
    network = ct.ReactorNet([reac])
    network.set_max_time_step(np.max(np.diff(v_of_t[0, :])))

    # Run the simulation
    time = 0
    states = ct.SolutionArray(gas, extra=['t'])
    states.append(reac.thermo.state, t=time)
    while time < np.min([end_time, v_of_t[0, -1]]):
        time = network.step()
        states.append(reac.thermo.state, t=time)

    # Get results
    times = states.t
    pressures = states.P
    temps = states.T
    target_concs = np.zeros((len(target_spcs), len(times)))
    for idx, target_spc in enumerate(target_spcs):
        target_concs[idx, :] = states.X[:, gas.species_index(target_spc)]
    end_gas = gas
    
    if run_rop: 
        pass
    else:
        rop = []

    return target_concs, pressures, temps, times, end_gas, rop


def rcm(temp, pressure, mix, gas, target_spcs, end_time, v_of_t):
    """ Runs a rapid compression machine simulation, with compression and heat
        loss incorporated via an input V(t) profile

        :param temp: reactor inlet temperature (Kelvin)
        :type temp: float
        :param pressure: reactor constant pressure (atm)
        :type pressure: float
        :param mix: initial species concentrations at reactor inlet
        :type mix: dct {spc1: conc1, spc2: ...} (example:
            {'H2': 0.5, 'O2': 0.5})
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param target_spcs: desired species concentrations
        :type target_spcs: list [spc1, spc2, ...]
        :param end_time: end time for the simulation (s); simulation will end
            at the first of either this time or the max time in v_of_t
        :type end_time: float
        :param v_of_t: specified volume (any units) vs time (s) profile
            to account for both piston compression and heat loss
        :type v_of_t: Numpy array of shape (2, num_times), where
            index 0 is time and index 1 is volume
        :return pressures: calculated pressure (Pa)
        :rtype: Numpy array of shape (num_times)
        :return times: time array corresponding to pressures (s)
        :rtype: Numpy array of shape (num_times)
    """

    def _piston_velocity(time):
        """ Gets the piston velocity at some time based on v_of_t

            :param time: individual time value (s)
            :type time: float
            :return vel: velocity of the piston (m/s)
            :rtype: float
        """
        dvdt = np.gradient(v_of_t[1, :], v_of_t[0, :])
        dxdt = dvdt / piston.area  # only for clarity (since A = 1 m^2)
        vel = np.interp(time, v_of_t[0, :], -1 * dxdt)

        return vel

    pressure = pressure * 101325  # atm to Pa
    gas.TPX = temp, pressure, mix

    # Normalize the volume profile by the first value
    v_of_t[1, :] = v_of_t[1, :] / v_of_t[1, 0]

    # Set up reactor, piston, and network
    reac = ct.IdealGasReactor(gas)
    reac.volume = 1  # m^3
    env = ct.Reservoir(gas)
    piston = ct.Wall(env, reac)
    piston.area = 1  # m^2 (this is the default, but put here for clarity)
    piston.heat_transfer_coeff = 0  # instead, considering heat trans. w/volume
    piston.set_velocity(_piston_velocity)
    network = ct.ReactorNet([reac])
    network.set_max_time_step(np.max(np.diff(v_of_t[0, :])))

    # Run the simulation
    time = 0
    states = ct.SolutionArray(gas, extra=['t'])
    states.append(reac.thermo.state, t=time)
    while time < np.min([end_time, v_of_t[0, -1]]):
        time = network.step()
        states.append(reac.thermo.state, t=time)

    # Get results
    times = states.t
    pressures = states.P
    target_concs = np.zeros((len(target_spcs), len(times)))
    for idx, target_spc in enumerate(target_spcs):
        target_concs[idx, :] = states.X[:, gas.species_index(target_spc)]
    end_gas = gas

    return target_concs, end_gas, rop, times, pressures


def pfr(temp, pressure, mix, gas, target_spcs, mdot, area, length,
        n_steps=2000):
    """ Runs a plug flow reactor simulation

        :param temp: reactor inlet temperature (Kelvin)
        :type temp: float
        :param pressure: reactor constant pressure (atm)
        :type pressure: float
        :param mix: initial species concentrations at reactor inlet
        :type mix: dct {spc1: conc1, spc2: ...}
            (example: {'H2': 0.5, 'O2': 0.5})
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param target_spcs: desired species concentrations
        :type target_spcs: list [spc1, spc2, ...]
        :param mdot: reactor mass flow rate (kg/s)
        :type mdot: float
        :param area: reactor cross-sectional area (m^2)
        :type area: float
        :param length: reactor length (m)
        :type length: float
        :param n_steps: approx. number of steps to take; used to guess timestep
            (note: will almost NEVER be the actual number of timesteps taken)
        :type n_steps: int
        :return target_concs: target_spcs mole fractions at each time/position
        :rtype: Numpy array of shape (num_targets, num_times)
        :return times: time values at each solver step
        :rtype: Numpy array of shape (num_times,)
        :return positions: position values at each solver step
        :rtype: Numpy array of shape (num_times,)
    """
    pressure = pressure * 101325  # atm to Pa
    gas.TPX = temp, pressure, mix

    # Create reactor and reactor network
    reac = ct.IdealGasConstPressureReactor(gas)
    network = ct.ReactorNet([reac])

    # Approximate a time step
    density, _ = gas.DP  # kg/m^3
    inlet_velocity = mdot / (density * area)  # m/s
    dt = (length / inlet_velocity) / n_steps  # s

    # Initialize time, position, velocity, and thermo states
    times = (np.arange(n_steps * 2) + 1) * dt  # longer than needed to be safe
    positions = np.zeros_like(times)
    velocities = np.zeros_like(times)
    states = ct.SolutionArray(reac.thermo)

    # Loop over each timestep
    end_idx = len(times)  # defining as backup; should be overwritten
    for idx, time in enumerate(times):
        network.advance(time)  # perform time integration
        velocities[idx] = mdot / area / reac.thermo.density
        positions[idx] = positions[idx - 1] + velocities[idx] * dt  # transform
        states.append(reac.thermo.state)

        # If the reactor length has been exceeded, exit the for loop
        if positions[idx] >= length:
            end_idx = idx
            break

    # Removed unused entries in time and position
    times = times[:(end_idx + 1)]
    positions = positions[:(end_idx + 1)]

    # Get results for target species
    target_concs = np.zeros((len(target_spcs), end_idx + 1))
    for idx, target_spc in enumerate(target_spcs):
        target_concs[idx, :] = states.X[:, gas.species_index(target_spc)]

    return target_concs, times, positions


def jsr(temp, pressure, mix, gas, target_spcs, res_time, vol, prev_concs=None,
        mdot=None, max_iter=30000):
    """ Runs a jet-stirred reactor simulation

        :param temp: reactor inlet temperature (Kelvin)
        :type temp: float
        :param pressure: reactor constant pressure (atm)
        :type pressure: float
        :param mix: initial species concentrations at reactor inlet
        :type mix: dct {spc1: conc1, spc2: ...} (example:
            {'H2': 0.5, 'O2': 0.5})
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param target_spcs: desired species concentrations
        :type target_spcs: list [spc1, spc2, ...]
        :param res_time: reactor residence time (s)
        :type res_time: float
        :param vol: volume of the reactor (m^3)
        :type vol: float
        :param prev_concs: species concentrations from previous solution at
            a similar condition
        :type prev_concs: numpy array of shape (num_spcs_in_mechanism,)
        :param max_iter: max number of timesteps to achieve steady state
        :type max_iter: int
        :return target_concs: outlet concentrations of target_spcs
        :rtype: Numpy array of shape (num_targets,)
        :return all_concs: outlet concentrations of all species in mech
        :rtype: Numpy array of shape (num_spcs_in_mechanism,)
    """

    assert res_time is None or mdot is None, (
        'Both res_time and mdot were provided; only one should be given.')

    pressure = pressure * 101325  # atm to Pa
    gas.TPX = temp, pressure, mix

    # Create reservoirs
    inlet = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)

    # Create reactor
    prev_concs_input = True
    if prev_concs is None:
        prev_concs_input = False
        prev_concs = mix
    gas.TPX = temp, pressure, prev_concs  # prev_concs to speed up convergence
    reac = ct.IdealGasReactor(gas, energy='off', volume=vol)

    # Set up devices
    pressure_valve_coeff = 0.01  # "conductance" of the pressure valve
    p_reg = ct.Valve(upstream=reac, downstream=exhaust, K=pressure_valve_coeff)
    # The mass flow controller condition depends on the inputs
    if res_time is not None:
        mfc = ct.MassFlowController(upstream=inlet, downstream=reac,
                                    mdot=reac.mass / res_time)
    elif mdot is not None:
        mfc = ct.MassFlowController(upstream=inlet, downstream=reac, mdot=mdot)
    else:  # if both res_time and mdot are None
        raise NotImplementedError('Either mdot or res_time must be given.')
    # print('res_time calc:\n', reac.mass/mdot)

    # Create reactor network (only the JSR in this case) and advance it to SS
    reac_net = ct.ReactorNet([reac])
    failure = False
    try:
        reac_net.advance_to_steady_state(max_steps=max_iter)
        all_concs = reac.thermo.X  # store output concentrations
    except ct._cantera.CanteraError:
        failure = True
        print(f'The solver failed at {temp} K for mechanism {gas.name}')
        # If no initial guess, next iteration will use initial mix as guess
        if prev_concs_input is False:
            all_concs = None
        # If an initial guess was input, return it for use in the next iteration
        else:
            all_concs = prev_concs

    # Get results for target species
    target_concs = np.zeros(len(target_spcs))
    for idx, target_spc in enumerate(target_spcs):
        if failure:
            target_concs[idx] = None
        else:
            if target_spc in gas.species_names:
                target_concs[idx] = all_concs[gas.species_index(target_spc)]
            else:  # if the target_spc isn't in the mechanism
                target_concs[idx] = None
    end_gas = gas

    return target_concs, all_concs, end_gas


def const_t_p(temp, pressure, mix, gas, target_spcs, end_time):
    """ Runs a constant-temperature, constant-pressure 0-D simulation

        :param temp: reactor inlet temperature (Kelvin)
        :type temp: float
        :param pressure: reactor constant pressure (atm)
        :type pressure: float
        :param mix: initial species concentrations at reactor inlet
        :type mix: dct {spc1: conc1, spc2: ...} (example:
            {'H2': 0.5, 'O2': 0.5})
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param target_spcs: desired species concentrations
        :type target_spcs: list [spc1, spc2, ...]
        :param end_time: end time for the simulation (s)
        :type end_time: float
        :return target_concs: target_spcs mole fractions at each time
        :rtype: Numpy array of shape (num_targets, num_times)
        :return times: array of time values (s)
        :rtype: Numpy array of shape (num_times,)
    """

    pressure = pressure * 101325  # atm to Pa
    gas.TPX = temp, pressure, mix

    # Setting energy to 'off' holds T constant
    reac = ct.IdealGasConstPressureReactor(gas, energy='off')
    network = ct.ReactorNet([reac])

    # Run the simulation
    time = 0
    states = ct.SolutionArray(gas, extra=['t'])
    states.append(reac.thermo.state, t=time)
    while time < end_time:
        time = network.step()
        states.append(reac.thermo.state, t=time)

    # Get results
    times = states.t
    num_targets = len(target_spcs)  # overly verbose for clarity
    num_times = len(times)
    target_concs = np.zeros(num_targets, num_times)
    for idx, target_spc in enumerate(target_spcs):
        target_concs[idx, :] = states.X[:, gas.species_index(target_spc)]

    return target_concs, times
