# Year-round Operation Parameter Table

Definitions for symbols used across the full report Year-round Operation section (Sections 18.1 through 18.6), including explicit aliases such as $\delta^h_{m,q}$ and $\delta^f_{m,q}$.

| Symbol | Definition |
|---|---|
| $\alpha$ | wire temperature coefficient of resistivity. |
| $\bar T_m$ | monthly Connecticut mean ambient from `climate\_data.monthlyMean\_C`. |
| $\bar T_{w,g}(d)$ | representative wire temperature for group $g$, day $d$. |
| $\beta_{f,m},\rho_{f,m},U_{f,m}^{\mathrm{cnt}}$ | freeze-count floor, freeze-count fractional remainder, and freeze-count uniform draw for month $m$. |
| $\beta_{h,m},\rho_{h,m},U_{h,m}^{\mathrm{cnt}}$ | heat-wave-day count floor, fractional remainder, and uniform draw for month $m$. |
| $\beta_{s,m},\rho_{s,m},U_{s,m}^{\mathrm{cnt}}$ | spell-start count floor, fractional remainder, and uniform draw for month $m$. |
| $\chi_{\mathrm{derive}}$ | derivation toggle from `deriveTargetFromRatedSpecs` ($1$ on, $0$ off). |
| $\Delta p_{\mathrm{avail},x}(\dot V)$ | available blower pressure at trial flow $\dot V$ from the linear fan-curve relation. |
| $\Delta p_{\mathrm{dist,common}}$ | heating common distribution losses (header and header-to-splitter path). |
| $\Delta p_{\mathrm{dist,private},g}$ | heating group-$g$ private distribution losses (splitter body, connector friction, connector-to-branch area change). |
| $\Delta p_{\mathrm{dist},c}$ | cooling distribution loss sum $\Delta p_{\mathrm{header}}+\Delta p_{\mathrm{header\rightarrow splitter}}+\Delta p_{\mathrm{splitter\,body}}+\Delta p_{\mathrm{connector\,fric}}+\Delta p_{\mathrm{connector\rightarrow branch}}$. |
| $\Delta p_{\mathrm{Ergun},s}$ | porous-bed backpressure increment over path length $L_s$ in segment $s$. |
| $\Delta p_{\mathrm{heater,tube},g},\Delta p_{\mathrm{heater,coil},g},\Delta p_{\mathrm{heater\rightarrow aer},g}$ | heater-tube friction, heater-coil obstruction, and heater-to-aeration area-change losses for group $g$. |
| $\Delta p_{\mathrm{op},x}(d)$ | operating pressure used for power at solved day-$d$ flow, bounded by system-required and available pressures. |
| $\Delta p_{\mathrm{path},g}$ | full group-$g$ private path pressure requirement in heating. |
| $\Delta p_{\mathrm{rated},x}$ | static pressure \emph{rise} reported at the same rated operating point as $\dot V_{\mathrm{rated},x}$; in code this comes from `ratedPressure\_inH2O` and is converted to Pa. |
| $\Delta p_{\mathrm{shutoff},x}$ | blower $x$ shutoff pressure, i.e., the maximum static pressure \emph{rise} at zero flow from the blower datasheet; in code this comes from `shutoffPressure\_inH2O` and is converted to Pa. |
| $\Delta p_{\mathrm{sys},x}(\dot V)$ | system-required pressure drop at trial flow $\dot V$ for blower path $x$. In cooling, $\Delta p_{\mathrm{sys},\mathrm{ab}}=\Delta p_{\mathrm{dist},c}+p_{\mathrm{in},c}$. In heating, $\Delta p_{\mathrm{sys},\mathrm{hb}}=\Delta p_{\mathrm{dist,common}}+\Delta p_{\mathrm{path,ctrl}}$. |
| $\Delta p_{\mathrm{tube},s}$ | in-tube Darcy friction drop across segment $s$. |
| $\Delta T_{f,m,q}$ | sampled anomaly for selected freeze event $q$ in month $m$ (bounded to $[\Delta T_{f,m}^{\mathrm{low}},\Delta T_{f,m}^{\mathrm{high}}]$). |
| $\Delta T_{f,m}$ | freeze anomaly below monthly mean |
| $\Delta T_{f,m}^{\mathrm{low}},\Delta T_{f,m}^{\mathrm{high}}$ | lower/upper freeze-anomaly bounds for month $m$. |
| $\Delta T_{h,m,q}$ | sampled anomaly assigned to spell $q$ in month $m$, bounded to $[\Delta T_{h,m}^{\mathrm{low}},\Delta T_{h,m}^{\mathrm{high}}]$. |
| $\Delta T_{h,m}$ | monthly mean spell anomaly above monthly mean |
| $\Delta T_{h,m}^{\mathrm{low}},\Delta T_{h,m}^{\mathrm{high}}$ | lower/upper spell-anomaly bounds (above $\bar T_m$) for month $m$, loaded from NOAA-derived monthly spell statistics. |
| $\delta^f_{m,q}$ | Alias for freeze-event anomaly magnitude; identical to $\Delta T_{f,m,q}$. |
| $\delta^h_{m,q}$ | Alias for heat-wave spell anomaly magnitude; identical to $\Delta T_{h,m,q}$. |
| $\delta_i(d)$ | circular day distance (Eq. 230) |
| $\dot m(d)$ | total moist-air mass flow through the spot-cooler branch, computed as $\rho_{\mathrm{air}}(T_\infty(d),p_\infty)\dot V_{\mathrm{ach},\mathrm{ab}}(d)$. |
| $\dot m(d),c_p(d)$ | day-$d$ cooling-branch mass flow and effective specific heat. |
| $\dot m_{\mathrm{dry}}(d)$ | dry-air mass flow, $\dot m(d)/(1+w_{\mathrm{in}}(d))$. |
| $\dot m_{\mathrm{evap},b,s}$ | evaporation mass flow rate in branch $b$, segment $s$. |
| $\dot m_{\mathrm{mean},b,s}$ | mean airflow through segment $s$ of branch $b$. |
| $\dot m_{\mathrm{rated}}(d)$ | rated mass flow at day-$d$ ambient properties. |
| $\dot m_{\mathrm{rel},b,s}$ | mass flow released through perforations into bed in branch $b$, segment $s$. |
| $\dot m_{\mathrm{remaining}}(p_{\mathrm{in}})$ | remaining mass-flow mismatch used in $R_{\mathrm{inlet}}$. |
| $\dot Q_{\mathrm{rated,sens}}(d)$ | rated sensible capacity used by the target-derivation path. |
| $\dot Q_{\mathrm{sensible}}(d)$ | sensible (dry-bulb) cooling rate, defined as $\dot m(d)c_p(d)\max(T_\infty(d)-T_{\mathrm{supply,actual}}(d),0)$, equal to $\dot Q_{\mathrm{spot,act}}(d)$ in this model. |
| $\dot Q_{\mathrm{spot,act}}(d)$ | actual spot-cooler sensible load after capacity cap. |
| $\dot Q_{\mathrm{spot,max}}(d)$ | available spot-cooler sensible capacity at day $d$. |
| $\dot Q_{\mathrm{spot,req}}(d)$ | requested sensible spot-cooling load at day $d$. |
| $\dot V_{\mathrm{ach},\mathrm{ab}}(d)$ | achieved cooling assist-blower flow from the blower/network residual solve. |
| $\dot V_{\mathrm{ach},x}(d)$ | achieved flow for blower $x$, equal to requested flow under disabled/valid-request-feasible conditions, otherwise $\dot V_{\mathrm{feas},x}(d)$. |
| $\dot V_{\mathrm{feas},x}(d)$ | helper feasible-flow variable, equal to $\max\mathcal{F}_x(d)$ when feasible set is nonempty, else $0$. |
| $\dot V_{\mathrm{free},x}$ | blower $x$ free-air flow (fan-curve intercept at zero pressure). |
| $\dot V_{\mathrm{rated},x}$ | blower $x$ rated volumetric flow from the manufacturer-rated operating point; in code this comes from `ratedFlow\_CFM` and is converted to SI units. |
| $\dot V_{\mathrm{rated}}$ | rated airflow converted to $\mathrm{m^3/s}$. |
| $\dot V_{\mathrm{req},x}(d)$ | requested blower volumetric flow for day $d$, blower $x$. |
| $\ell_{h,m,q}$ | Heat-wave spell length (days) for spell $q$ in month $m$. |
| $\eta_x$ | total blower efficiency used in $P_{\mathrm{shaft}}/\eta$ conversion for blower $x$. |
| $\hat y(d)$ | kernel-regressed daily estimate generated from one selected $y_i$ pass at day $d$. |
| $\mathbf{A}(d),\mathbf{T}(d),\mathbf{b}(d)$ | day-$d$ two-node linear-system matrix, unknown temperature vector, and right-hand-side vector. |
| $\mathcal{B}_{h,m,q}$ | Set of day-of-month indices occupied by heat-wave spell $q$ in month $m$. |
| $\mathcal{B}_{h,m}^{\mathrm{all}}$ | Union of all heat-wave day blocks in month $m$: $\bigcup_q \mathcal{B}_{h,m,q}$. |
| $\mathcal{F}_x(d)$ | feasible-flow set for blower $x$ on day $d$, containing trial flows in $[0,\dot V_{\mathrm{req},x}(d)]$ with $R_x(\dot V)\le 0$. |
| $\mathcal{I}_{f,m}^{(k)}$ | freeze-day candidates remaining at step $k$, defined by $\mathcal{I}_{f,m}^{(k)}=\{1,\dots,L_m\}\setminus\mathcal{S}_{f,m}^{(k-1)}$. |
| $\mathcal{I}_{h,m}^{(k)}$ | heat-wave-start candidates remaining at step $k$, defined by $\mathcal{I}_{h,m}^{(k)}=\{1,\dots,L_m\}\setminus\mathcal{S}_{h,m}^{(k-1)}$. |
| $\mathcal{S}_{f,m}$ | set of freeze-day indices already selected in month $m$. |
| $\mathcal{S}_{f,m}^{(k)}$ | freeze-day selected set after $k$ recursive selections ($\mathcal{S}_{f,m}^{(0)}=\varnothing$). |
| $\mathcal{S}_{h,m}$ | set of heat-wave start-day indices selected in month $m$. |
| $\mathcal{S}_{h,m}^{(k)}$ | heat-wave-start selected set after $k$ recursive selections ($\mathcal{S}_{h,m}^{(0)}=\varnothing$). |
| $\mathrm{COP}$ | coefficient of performance used to convert sensible load to electric power. |
| $\mu_{f,m}$ | mean freeze day-of-month |
| $\mu_{h,m}$ | mean heat-wave start day |
| $\rho_{20}$ | wire resistivity at $20\,^\circ\mathrm{C}$. |
| $\rho_{\mathrm{air}}(T,p)$ | moist-air density function used to map rated airflow to rated mass flow. |
| $\sigma(d)$ | regressed daily spread from the same kernel pass on $\sigma_m$. |
| $\sigma_m$ | monthly ambient spread from `climate\_data.monthlyStd\_C`. |
| $\sigma_{\Delta T,f,m}$ | freeze-anomaly spread parameter for month $m$. |
| $\sigma_{\Delta T,h,m}$ | configured or fallback standard deviation used for spell-anomaly sampling in month $m$. |
| $\varepsilon,d_p,L_s,A_{\mathrm{trib},s}$ | bed porosity, representative particle diameter, escape-path length, and tributary porous area in segment $s$. |
| $\varnothing$ | Empty set. Also used with set operators: $A\setminus B$ (difference), $A\cup\{x\}$ (union/add), and $\|A\|$ (set cardinality). |
| $A_{f,m}(r)$ | day-of-month freeze-anomaly field; $A_{f,m}(r)=\Delta T_{f,m,q}$ when day $r$ is the $q$-th selected freeze day, and $A_{f,m}(r)=0$ otherwise. |
| $A_{h,m}(r)$ | day-of-month spell-anomaly field for month $m$, obtained by allocating each spell magnitude $M_{h,m,q}$ to all days in spell block $\mathcal{B}_{h,m,q}$. |
| $b$ | aeration branch index ($1,\dots,N_b$). |
| $C(d)$ | daily cost (Eq. 250) |
| $c_{\mathrm{elec}}$ | electricity price (Eq. 250) |
| $c_{p,b,s}$ | specific heat of segment airflow in branch $b$, segment $s$. |
| $d$ | day-of-year evaluation index ($d=1,\dots,P$). |
| $D(d)$ | Active duty fraction for day $d$: equals $D_h(d)$ on heating days or $D_c(d)$ on cooling days. |
| $D_h(d),D_c(d)$ | mode-specific clipped duty fractions. |
| $D_h^*(d),D_c^*(d)$ | required duty fraction (Eq. 243) |
| $d_i$ | day-of-year reference day associated with reference point $i$ (month midpoint in this implementation). |
| $d_m$ | month-midpoint day-of-year for month $m$. With one reference point per month, $d_m$ and $d_i$ refer to the same reference-day map. |
| $d_{f,m,k}$ | freeze day selected at recursion step $k$. |
| $d_{h,m,k}$ | heat-wave start day selected at recursion step $k$. |
| $d_{h,m,q}^{\mathrm{raw}}$ | Raw ordered start day for heat-wave spell $q$ before feasibility projection. |
| $E_h(d),E_c(d)$ | daily mode-specific kWh (Eq. 248) |
| $E_{\mathrm{day}}$ | total daily kWh (Eq. 249) |
| $e_{h,m,q}$ | Earliest feasible start day for heat-wave spell $q$. |
| $F_m$ | expected freeze days/month |
| $g$ | heater-tube group index; $n_{\mathrm{tube},g}$ is the number of identical tubes in group $g$. |
| $g_{\mathrm{ctrl}}$ | controlling heating group (the group with the largest $\Delta p_{\mathrm{path},g}$). |
| $g_{\mathrm{std}}(z),G_{\mathrm{std}}(z)$ | standard-normal density and CDF used in the explicit truncated-normal equations. |
| $h$ | Gaussian kernel bandwidth in days. |
| $H_m$ | expected heat-wave days/month |
| $h_{fg,b,s}$ | latent heat of vaporization used in branch $b$, segment $s$. |
| $I_{\mathrm{tube},g}(d)$ | current through one tube/string in group $g$. |
| $L_m$ | number of calendar days in month $m$. |
| $L_{w,g},A_{w,g}$ | wire length and wire cross-sectional area for group $g$. |
| $m$ | month index ($m=1,\dots,12$). |
| $m_{\mathrm{water}}(d)$ | duty-scaled daily water replacement (Eq. 251) |
| $M_{f,m}$ | freeze-anomaly monthly center value for month $m$ (mean freeze anomaly below $\bar T_m$). |
| $M_{h,m,q}$ | spell-magnitude variable for spell $q$ in month $m$; in code-equivalent notation $M_{h,m,q}\equiv \Delta T_{h,m,q}$. |
| $N_f,N_s,N_h$ | realized annual totals after month-wise stochastic rounding (freeze days, spell starts, and heat-wave days). |
| $n_{f,m}$ | allocated integer freeze days (Eq. 231) |
| $N_{s,g}$ | number of tube elements in series in group $g$. |
| $n_{s,m},n_{h,m}$ | allocated starts and days (Eq. 235) |
| $P$ | annual period length in days ($P=365$ in this implementation). |
| $p_\infty$ | day ambient absolute pressure used in air-property evaluation. |
| $P_h(T_\infty(d))$ | heating-mode electrical power from the day-ambient heating re-solve. |
| $P_{\mathrm{assist\,blower}}(d)$ | assist-blower electric power on day $d$ from the cooling re-solve. |
| $P_{\mathrm{hb},\mathrm{elec}}(d)$ | heating-blower electrical power at day $d$. |
| $P_{\mathrm{idle},x}$ | baseline low-flow electrical draw for blower $x$ (motor/controller/parasitics), from `idlePower\_W`; in implemented logic it is nonnegative and only applied when blower is enabled and $\dot V_{\mathrm{ach},x}>0$, otherwise blower electrical power is set to $0$. |
| $p_{\mathrm{in,aer},g}$ | group-$g$ aeration inlet gauge pressure from the same closed-end release/backpressure closure. |
| $p_{\mathrm{in},c}$ | cooling aeration-branch inlet gauge pressure from closed-end inlet solve. |
| $P_{\mathrm{in}}(d)$ | spot-cooler electrical input power used in COP definitions; in this model $P_{\mathrm{in}}(d)\equiv P_{\mathrm{spot\,cooler}}(d)$. |
| $p_{\mathrm{out,bnd}}$ | bed outlet-boundary gauge pressure (open-top/vent/offset closure). |
| $P_{\mathrm{rated},x}$ | motor rated electrical power cap for blower $x$. |
| $P_{\mathrm{shaft},x}(d)$ | shaft (fluid) power for blower $x$, $P_{\mathrm{shaft},x}=\dot V_{\mathrm{ach},x}\Delta p_{\mathrm{op},x}$. |
| $P_{\mathrm{spot\,cooler}}(d)$ | spot-cooler electric power on day $d$ from the cooling re-solve. |
| $P_{\mathrm{tube},g}(d)$ | electrical power of one tube/string in group $g$. |
| $P_{c,\mathrm{elec}}(d)$ | cooling electric power (Eq. 242) |
| $P_{h,\mathrm{coil}}(d)$ | total heater-coil electrical power from all tube groups at day $d$. |
| $P_{h,\mathrm{elec}}(d)$ | heating electric power (Eq. 241) |
| $P_{h,\mathrm{elec}}(d),P_{c,\mathrm{elec}}(d)$ | mode electric powers used in daily energy equations, where $P_{h,\mathrm{elec}}(d)=P_h(T_\infty(d))$ and $P_{c,\mathrm{elec}}(d)=P_{\mathrm{spot\,cooler}}(d)+P_{\mathrm{assist\,blower}}(d)$. |
| $P_{x,\mathrm{elec}}(d)$ | blower-$x$ electrical power at day $d$, computed by efficiency model (or flow-fraction fallback), then capped by $P_{\mathrm{rated},x}$. |
| $q$ | heat-wave spell index within month $m$, ordered by projected start day. |
| $Q_{\mathrm{bed}}(d)$ | applied bed load (Eq. 246) |
| $Q_{\mathrm{evap},b,s}$ | latent heat sink due to evaporation in branch $b$, segment $s$. |
| $Q_{\mathrm{jet},b,s}$ | sensible heat carried by released jet flow relative to bed reference temperature. |
| $Q_{\mathrm{req},c}(d)$ | required cooling load (Eq. 240) |
| $Q_{\mathrm{req},h}(d)$ | required heating load (Eq. 239) |
| $Q_{\mathrm{to\,bed},c}(T_\infty(d))$ | signed net heat-to-bed returned by the cooling reference-point re-solve at day ambient $T_\infty(d)$ (typically negative when cooling removes heat from the bed). |
| $Q_{\mathrm{to\,bed},h}(T_\infty(d))$ | signed net heat-to-bed returned by the heating reference-point re-solve at day ambient $T_\infty(d)$ (positive when the heater adds heat to the bed). |
| $Q_{\mathrm{wall},b,s}$ | sensible heat transfer from branch airflow to bed-side walls in branch $b$, segment $s$. |
| $Q_{c,\mathrm{cap}}(d)$ | cooling capacity at day ambient (Eq. 242) |
| $Q_{h,\mathrm{cap}}(d)$ | heating capacity at day ambient (Eq. 241) |
| $Q_{h,\mathrm{cap}}(d),Q_{c,\mathrm{cap}}(d)$ | positive heating and cooling capacities used in duty equations, with $Q_{h,\mathrm{cap}}(d)=Q_{\mathrm{to\,bed},h}(T_\infty(d))$ and $Q_{c,\mathrm{cap}}(d)=-Q_{\mathrm{to\,bed},c}(T_\infty(d))$. |
| $r$ | generic day-of-month index, $r\in\{1,\dots,L_m\}$. |
| $r_m(d)$ | day-of-month index corresponding to day-of-year $d$ in month $m$. |
| $R_x(\dot V)$ | blower/system residual for blower $x$, $R_x=\Delta p_{\mathrm{sys},x}-\Delta p_{\mathrm{avail},x}$. |
| $R_{\mathrm{bed\,ref}}(T_{\mathrm{ref}})$ | bed-reference closure residual in cooling branch. |
| $R_{\mathrm{inlet}}(p_{\mathrm{in}})$ | closed-end inlet mass-balance residual $\dot m_{\mathrm{rem},N_{\mathrm{seg}}+1}(p_{\mathrm{in}})$. |
| $R_{\mathrm{tube},g}(d)$ | electrical resistance of one tube in group $g$ at day-$d$ operating temperature. |
| $R_{f,m}^{(k)}(r)$ | deterministic ranking value for freeze-day candidate day $r$ at step $k$; larger $R$ means higher selection priority. |
| $R_{h,m}^{(k)}(r)$ | deterministic ranking value for heat-wave-start candidate day $r$ at step $k$; larger $R$ means higher selection priority. |
| $s$ | axial segment index ($1,\dots,N_{\mathrm{seg}}$). |
| $S_m$ | expected heat-wave spell starts/month |
| $s_{f,m},r_{f,m}$ | freeze-day Gaussian spread and repulsion length scales for month $m$. |
| $s_{h,m,q}$ | Feasible heat-wave spell start day for spell $q$ in month $m$ after non-overlap and month-bound clipping. |
| $s_{h,m},r_{h,m}$ | heat-wave start Gaussian spread and repulsion length scales for month $m$. |
| $T_\infty(d)$ | final day ambient passed to thermal load and capacity re-solves. |
| $T_b(d),T_t(d)$ | bottom/top bed-node temperatures from the two-node steady solve under $Q_{\mathrm{bed}}(d)$ and $T_\infty(d)$. |
| $T_b,T_t$ | bottom and top bed-node temperatures used by the closure. |
| $T_m(d)$ | month-local daily ambient profile used during event overlays inside month $m$. |
| $T_{\infty,\mathrm{base}}(d)$ | regressed daily ambient baseline from periodic Gaussian kernel regression on $\bar T_m$. |
| $T_{\infty,\mathrm{heatRaw}}(d)$ | optional diagnostic track exported by code; mathematically, Eq. 18.3-G is written directly in combined form for $T_\infty(d)$. |
| $T_{\infty,\mathrm{preH}}(d)$ | intermediate daily ambient after freeze-day cold overrides (input to 18.3 heat-wave overlay). |
| $T_{\mathrm{air,in},b,s},T_{\mathrm{air,out},b,s}$ | segment inlet/outlet air temperatures. |
| $T_{\mathrm{bed,ref}}$ | branch-level bed reference temperature used for local segment balances. |
| $T_{\mathrm{ref}}$ | trial bed reference temperature in the cooling closure. |
| $T_{\mathrm{rel},b,s}$ | release/jet air temperature to the bed in branch $b$, segment $s$. |
| $T_{\mathrm{supply,actual}}(d)$ | actual spot-cooler supply-air temperature after capacity capping. |
| $T_{\mathrm{target,cfg}}$ | configured spot-cooler target from `targetSupplyAir\_C`. |
| $T_{\mathrm{target,derived}}(d)$ | target supply temperature derived from rated sensible capacity and rated mass flow. |
| $T_{\mathrm{target}}(d)$ | effective day-$d$ spot-cooler target used in $\dot Q_{\mathrm{spot,req}}(d)$. |
| $T_{h,\mathrm{set}},T_{c,\mathrm{set}}$ | Resolved heating and cooling setpoints used by the year-round solver. |
| $u_s$ | superficial porous velocity used in Ergun pressure-drop term for segment $s$. |
| $u_{h,m,q}$ | Latest feasible start day for heat-wave spell $q$. |
| $V$ | applied electrical voltage across each series string. |
| $V_{\mathrm{rated,CFM}}$ | spot-cooler rated airflow from configuration. |
| $w_i(d)$ | Gaussian kernel weight applied to reference point $i$ when evaluating day $d$. |
| $w_{\mathrm{in}}(d)$ | inlet humidity ratio used to convert moist-air to dry-air basis. |
| $x\in\{\mathrm{hb},\mathrm{ab}\}$ | blower index ($\mathrm{hb}$: heating blower, $\mathrm{ab}$: cooling assist blower). |
| $y_i$ | Reference monthly value used in the active kernel pass; either $\bar T_i$ (baseline pass) or $\sigma_i$ (spread pass). |
