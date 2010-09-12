BARDIR=../bar

bar: $(BARDIR)/rphi2d.ps $(BARDIR)/rphi1d.ps $(BARDIR)/detecta.ps \
	$(BARDIR)/detectb.ps $(BARDIR)/rphi1d2.ps \
	$(BARDIR)/distuncertain.ps $(BARDIR)/df.ps \
	$(BARDIR)/barstrength.ps $(BARDIR)/slope.ps \
	$(BARDIR)/apogee.ps \ $(BARDIR)/vlosuncertain.ps \
	$(BARDIR)/rolr.ps

$(BARDIR)/rphi2d.ps: bar_figures.py calc_veldist_2d.py
	$(PYTHON) bar_figures.py $@

$(BARDIR)/rphi1d.ps: bar_figures.py calc_veldist_1d.py
	$(PYTHON) bar_figures.py $@ --vlos

$(BARDIR)/rphi1d2.ps: bar_figures.py calc_veldist_1d.py
	$(PYTHON) bar_figures.py $@ --vlos --otherside

$(BARDIR)/detecta.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@

$(BARDIR)/detectb.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@ --skipCenter=0.15

$(BARDIR)/barstrength.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --barstrength

$(BARDIR)/df.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --df

$(BARDIR)/rolr.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --rolr

$(BARDIR)/slope.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --slope

$(BARDIR)/distuncertain.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --convolve

$(BARDIR)/vlosuncertain.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --vrconvolve

$(BARDIR)/apogee.ps: bar_figures.py
	$(PYTHON) bar_figures.py $@ --apogee


