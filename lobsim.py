from sys import platform

from cffi import FFI
import plotly.graph_objects as go


def run(candleNum):
    ffi = FFI()
    ffi.cdef("""
        void sim(int32_t stf, int32_t simulation_time, int32_t candle_size_ms, double tick_size, int32_t k, double p_ref,
        int32_t, double theta, double theta_reinit, double * o, double * h, double * l, double * c ,int32_t candle_num,
        double * calib);""")

    if "linux" in platform:
        pre = "lib"
        ext = ".so"
    else:  # assumes Windows
        pre = ""
        ext = ".dll"
    lobsim = ffi.dlopen('target/release/' + pre + 'lobsim' + ext)
    tick_size = 0.00001
    k = 3
    stf = 1

    assert candleNum % 1 == 0
    candle_size_ms = 1000 * 60 * 5
    simulation_time = 1000 * 60 * 5 * candleNum
    resize = 1
    p_ref = 1.230495
    theta, theta_reinit = 0.3, 0.5
    cdat = "double[" + str(candleNum) + "]"
    o = ffi.new(cdat)
    h = ffi.new(cdat)
    l = ffi.new(cdat)
    c = ffi.new(cdat)
    calib = ffi.new("double[2]")
    lobsim.sim(stf, simulation_time, candle_size_ms, tick_size, k, p_ref, resize, theta, theta_reinit, o, h, l, c,
               candleNum, calib)
    print("simulation complete, plotting...")
    candlesticks = go.Candlestick(open=list(o), high=list(h), low=list(l), close=list(c))
    fig = go.Figure(data=[candlesticks])
    fig.show()


if __name__ == "__main__":
    run(100)
