import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from iapws import IAPWS97
import numpy as np

st.title('Курсовая работа ')
st.subheader('Конюхова П.O. ФПэ-01-19 Вариант 6 ')
st.caption('Ne = 218 МВт')
st.caption('p0 = 13.1 МПа')
st.caption('t0 = 548 C')
st.caption('ppp = 2.489 МПа')
st.caption('tpp = 550 C')
st.caption('pk = 4.8 кПа')
st.caption('z = 7 шт.')
tpv = st.slider('Диапазон значений tpv = ',  min_value = 230, max_value = 240, step = 1)
tpv = tpv + 0.01

Ne = 218e6
p0 = 13.1e6
t0 = 548
T0 = t0 + 273.15
ppp = 2.49e6
tpp = 550
Tpp = tpp+273.15
pk = 4.8e3
tpv = list(np.arange(230, tpv, 1))
Tpv = [t+273.15 for t in tpv]

d_p0 = 0.05
d_p_pp = 0.08
d_p = 0.03

delta_p_0 = p0*d_p0
delta_p_pp = ppp*d_p_pp

def Calculate_G0_Gk(N_e, p_0, T_0, p_pp, T_pp, p_k, T_pv):
    d_p0 = 0.05
    d_p_pp = 0.08
    d_p = 0.03
    point_0 = IAPWS97(P=p_0*10**(-6),T=T_0)
    s_0 = point_0.s
    h_0 = point_0.h
    v_0 = point_0.v

    p_0_ = p_0-0.05*p_0
    point_p_0_ = IAPWS97(P=p_0_*10**(-6),h=h_0)
    t_0_ = point_p_0_.T-273.15
    s_0_ = point_p_0_.s
    v_0_ = point_p_0_.v
    p_1t = p_pp+0.1*p_pp
    point_1t = IAPWS97(P=p_1t*10**(-6),s=s_0)
    t_1t = point_1t.T-273.15
    h_1t = point_1t.h
    v_1t = point_1t.v
    point_pp = IAPWS97(P=p_pp*10**(-6),T=T_pp)
    h_pp = point_pp.h
    s_pp = point_pp.s
    v_pp = point_pp.v
    H_0 = h_0-h_1t
    eta_oi = 0.85
    H_i_cvd = H_0*eta_oi
    h_1 = h_0 - H_i_cvd
    point_1 = IAPWS97(P = p_1t*10**(-6),h = h_1)
    s_1 = point_1.s
    T_1 = point_1.T
    v_1 = point_1.v
    p_pp_ = p_pp - 0.03*p_pp
    point_pp_ = IAPWS97(P=p_pp_*10**(-6),h = h_pp)
    s_pp_ = point_pp_.s
    v_pp_ = point_pp_.v
    point_kt = IAPWS97(P = p_k*10**(-6),s = s_pp)
    T_kt = point_kt.T
    h_kt = point_kt.h
    v_kt = point_kt.v
    s_kt = s_pp
    H_0_csdcnd = h_pp-h_kt
    eta_oi = 0.85
    H_i_csdcnd = H_0_csdcnd*eta_oi
    h_k = h_pp - H_i_csdcnd
    point_k = IAPWS97(P = p_k*10**(-6), h = h_k)
    T_k = point_k.T
    s_k = point_k.s
    v_k = point_k.v
    point_k_v = IAPWS97(P = p_k*10**(-6),x=0)
    h_k_v = point_k_v.h
    s_k_v = point_k_v.s
    eta_oiI = (h_1-h_0)/(h_1t-h_0)
    p_pv = 1.4*p_0
    point_pv = IAPWS97(P = p_pv*10**(-6),T=T_pv)
    h_pv = point_pv.h
    s_pv = point_pv.s
    ksi_pp_oo = 1-(1-(T_k*(s_pp-s_k_v))/((h_0-h_1t)+(h_pp-h_k_v)))/(1-(T_k*(s_pp-s_pv))/((h_0-h_1t)+(h_pp-h_pv)))
    T_0_ = 374.2+273.15
    T_ = (point_pv.T - point_k.T) / (T_0_ - point_k.T)
    if T_ <= 0.636363636:
        ksi1 = -1.53 * T_ ** 2 + 2.1894 * T_ + 0.0048
    elif 0.636363636 < T_ <= 0.736363636:
        ksi1 = -1.3855 * T_ ** 2 + 2.0774 * T_ + 0.0321
    elif 0.736363636 < T_ <= 0.863636364:
        ksi1 = -2.6535 * T_ ** 2 + 4.2556 * T_ - 0.8569
    else:
        ksi1 = 21321
    ksi = ksi1
    ksi_r_pp = ksi*ksi_pp_oo
    eta_ir = (H_i_cvd+H_i_csdcnd)/(H_i_cvd+(h_pp-h_k_v))*1/(1-ksi_r_pp)
    H_i = eta_ir*((h_0-h_pv)+(h_pp-h_1))
    eta_m = 0.994
    eta_eg = 0.99
    G_0 = N_e/(H_i*eta_m*eta_eg*(10**3))
    G_k = N_e/((h_k-h_k_v)*eta_m*eta_eg*(10**3))*(1/eta_ir-1)
    return eta_ir, G_0, G_k

eta_ir_, G_0_, G_k_=[],[],[]
for T in Tpv:
    x = Calculate_G0_Gk(N_e=Ne, p_0=p0, T_0=T0, p_pp=ppp, T_pp=Tpp, p_k=pk, T_pv=T)
    eta_ir_.append(x[0])
    G_0_.append(x[1])
    G_k_.append(x[2])

 eta_ir_1=np.around(eta_ir_, decimals=3)
 G_0_1=np.around(G_0_, decimals=3)
 G_k_1=np.around(G_k_, decimals=3)


st.subheader(" КПД ")
st.text('Зависимость КПД от температуры Tpv')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot( [t for t in tpv ],[ Calculate_G0_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = pk, T_pv = t) for t
in Tpv ], )
ax.scatter( [t for t in tpv ],[ Calculate_G0_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = pk, T_pv = t) for t in Tpv ])
plt.xlabel("Tpv, C")
plt.ylabel("η, КПД")
st.pyplot(fig)


st.subheader(" Таблица полученных данных ")
st.text(' Значения КПД и расходов')
kpd = {'Температура, C': [ t for t in tpv],
'КПД, %': [100*eta_ir_1],
'Расход пара на входе в турбину, кг/с': [Calculate_G0(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = pk, T_pv = t) for t in Tpv ],
'Расход пара на входе в конденсатор, кг/c': [Calculate_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = pk, T_pv = t) for t in Tpv ],
}
df = pd.DataFrame(data=kpd)
st.table(df)

st.subheader(" График процесса расширения пара в турбине ")

fig = plt.figure()
ax = fig.add_subplot(111)
point_0 = IAPWS97(P=p0*1e-6, T=T0)
p_0_d = p0 - delta_p_0
point_0_d = IAPWS97(P=p_0_d*1e-6, h=point_0.h)
p_1t = ppp + delta_p_pp
point_1t = IAPWS97(P=p_1t*10**(-6), s=point_0.s)
H_01 = point_0.h - point_1t.h
kpd_oi = 0.85
H_i_cvd = H_01 * kpd_oi
h_1 = point_0.h - H_i_cvd
point_1 = IAPWS97(P=p_1t*1e-6, h=h_1)
point_pp = IAPWS97(P=ppp*1e-6, T=Tpp)
p_pp_d = ppp - delta_p_pp
point_pp_d = IAPWS97(P=p_pp_d*1e-6, h=point_pp.h)
point_kt = IAPWS97(P=pk*1e-6, s=point_pp.s)
H_02 = point_pp.h - point_kt.h
kpd_oi = 0.85
H_i_csd_cnd = H_02 * kpd_oi
h_k = point_pp.h - H_i_csd_cnd
point_k = IAPWS97(P=pk*1e-6, h=h_k)

s_0 = [point_0.s-0.05,point_0.s,point_0.s+0.05]
h_0 = [IAPWS97(P = p0*1e-6,s = s_).h for s_ in s_0]
s_1 = [point_0.s-0.05,point_0.s,point_0.s+0.18]
h_1 = [IAPWS97(P=p_1t*1e-6, s = s_).h for s_ in s_1]
s_0_d = [point_0_d.s-0.05, point_0_d.s, point_0_d.s+0.05]
h_0_d = h_0
s_pp = [point_pp.s-0.05,point_pp.s,point_pp.s+0.05]
h_pp = [IAPWS97(P=ppp*1e-6, s=s_).h for s_ in s_pp]
s_k = [point_pp.s-0.05,point_pp.s,point_pp.s+0.8]
h_k = [IAPWS97(P=pk*1e-6, s=s_).h for s_ in s_k]
s_pp_d = [point_pp_d.s-0.05,point_pp_d.s,point_pp_d.s+0.05]
h_pp_d = h_pp

plt.plot([point_0.s,point_0.s,point_0_d.s,point_1.s],[point_1t.h,point_0.h,point_0.h,point_1.h],'-ob')
plt.plot([point_pp.s,point_pp.s,point_pp_d.s,point_k.s],[point_kt.h,point_pp.h,point_pp.h,point_k.h],'-om')
plt.plot(s_0,h_0, '-g')
plt.plot(s_1,h_1, '-g')
plt.plot(s_0_d,h_0_d, '-g')
plt.plot(s_pp,h_pp, '-g')
plt.plot(s_k,h_k, '-g')
plt.plot(s_pp_d,h_pp_d, '-g')

for x, y, ind in zip([point_pp.s, point_k.s], [point_pp.h, point_k.h], ['{пп}', '{к}']):
   plt.text(x-0.35, y+40, '$h_' + ind + ' = %.2f $'%y)
for x, y, ind in zip([point_kt.s, point_pp_d.s], [point_kt.h, point_pp_d.h], ['{кт}', '{ппд}']):
   plt.text(x+0.03, y+40, '$h_' + ind + ' = %.2f $'%y)

for x, y, ind in zip ([point_0.s, point_1.s], [point_0.h, point_1.h], ['{0}', '{1}']):
   plt.text(x-0.1, y+100, '$h_' + ind + ' = %.2f $'%y)

for x, y, ind in zip([point_1t.s, point_0_d.s], [point_1t.h, point_0_d.h], ['{1т}', '{0д}']):
   plt.text(x+0.03, y-60, '$h_' + ind + ' = %.2f $'%y)
plt.title("h - s диаграмма")
plt.xlabel("s, кДж/(кг*С)")
plt.ylabel("h, кДж/кг")
plt.grid(True)
st.pyplot(fig)
