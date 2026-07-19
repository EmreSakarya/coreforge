# CoreForge — Teori ve Yöntem

Kod sisteminin çözdüğü denklemler, sayısal şemalar ve model
kalibrasyonlarının tam dökümü.

## 1 · Nötronik: çok-gruplu difüzyon özdeğer problemi

$$-\nabla\cdot(D_g\nabla\phi_g) + \Sigma_{r,g}\,\phi_g
 = \frac{\chi_g}{k_{eff}}\sum_{g'}\nu\Sigma_{f,g'}\phi_{g'}
 + \sum_{g'\neq g}\Sigma_{s,g'\to g}\,\phi_{g'},\qquad
 \Sigma_{r,g}=\Sigma_{a,g}+\sum_{g'\neq g}\Sigma_{s,g\to g'}$$

- Grup sayısı serbest (1–20, arayüzde 8), saçılma matrisi tam
  (up-scatter dâhil).
- **Ayrıklaştırma:** hücre-merkezli sonlu farklar; 2-D'de 5-nokta,
  3-D'de 7-nokta şablonu. Ara yüzey difüzyon katsayısı **harmonik
  ortalama**: c = 2·D₁D₂/(D₁+D₂)/Δ².
- **Sınır koşulları** (her yüz bağımsız): yansıtıcı (J=0) veya **Robin
  vakum** J_out = γ·φ_yüzey; FD biçimi c_vac = 2γD/[Δ(γΔ+2D)].
  γ=0.4692 (transport-düzeltmeli ekstrapolasyon; 0.5 = Marshak).
- **Özdeğer:** güç (fisyon kaynağı) iterasyonu; k, kaynak oranından
  güncellenir. Yakınsama çifte kriter: |Δk|<10⁻⁷ VE noktasal kaynak
  değişimi <10⁻⁵.
- **İç çözüm:** grup başına NINNER (=4) **red-black SOR** süpürmesi
  (ω=1.6); renk = (i+j+k) paritesi → aynı renk hücreleri bağımsız →
  **OpenMP** ile satır düzeyinde paralel. Tüm 7-nokta katsayıları
  iterasyon öncesi bir kez hesaplanır.
- **3-D temsil:** malzeme haritası katman aralıklarına atanır
  (`MAP l1 l2`); kısmi batmış kontrol çubukları ve eksenel
  reflektörler tam geometriyle modellenir. NZ=1 durumu 2-D motorla
  bit-düzeyinde özdeştir (verify: 2D ≡ extrude-3D, fark 0.01 pcm).

## 2 · Yakıt tasarımcısı (xslib): fizikselden makroskopiğe

Sayı yoğunlukları gerçek yoğunluk/atom kütlelerinden; **termal grup
gerçek 2200 m/s mikroskopik verileriyle** (U-235: σf 582.6 b, σa
680.9 b, ν 2.432; U-238 2.68 b; H 0.3326 b; doğal B 759 b; Pu-239
1017.9/747.4 b; Pu-241 1374.9/1012.3 b; Xe-135 2.65 Mb; Sm-149
40.14 kb) tek bir spektrum çarpanı S_TH altında — böylece nüklitler
arası GÖRELİ fizik (ör. boron/yakıt rekabeti) korunur. Hızlı grup,
yavaşlama (Σ₁→₂ = C_SD·N_H) ve transport sabitleri **bir kez** nominal
taze PWR hücresine (e=3.1 w/o, 0 ppm → D=[1.43,0.38], Σ₁₂=0.0165,
k∞=1.320) kalibre edilir ve bir daha değiştirilmez; zenginlik/boron/
yoğunluk cevabı fizikten gelir. Toplu fisyon ürünü çifti (90 b termal /
4 b hızlı), 0-D reaktivite-sınırlı yanmayı literatür bandına (%3.1 için
~40 MWd/kgU) oturtacak şekilde kalibredir — gerekçesi kod içinde
belgelidir.

**Ters çözücü:** hedef makroskopik setin νΣf₂'sinden zenginlik,
Σa₂'sinden boron sıralı bisection ile bulunur (her ikisi kanıtlanmış
monoton); kalan sabitlerin farkları rezidü olarak raporlanır.

## 3 · Tükenim (burnup): quasi-statik çevrim

Adım döngüsü: akı çöz → verilen özgül güce [W/gU] mutlak normalize et →
blok başına Bateman zincirini entegre et → kesitleri N'den yeniden kur.

Zincir: U-235 (yutulma), U-238 →(yakalama)→ Pu-239 →(yakalama)→
Pu-240 →(yakalama)→ Pu-241; fisyondan **denge Xe-135**
(N_Xe = γF/(λ+σφ)) ve **denge Sm-149**; toplu FP çifti birikimli.
Entegrasyon: adım içinde donmuş reaksiyon hızlarıyla RK4 (40 alt-adım).
Enerji muhasebesi: BU [MWd/kgU] = SP·t/1000; doğrulama — istenen çevrim
yanması ile yanma haritası ortalaması birebir örtüşür.

**Letdown / otomatik çevrim sonu:** her adımda kritik boron bisection
ile bulunur; boron 0 ppm'e inip k<1 olduğunda çevrim
"reactivity-limited" olarak biter → ulaşılabilir yanma + EFPD raporlanır.

## 4 · Transient: nokta reaktör kinetiği

$$\frac{dn}{dt}=\frac{\rho(t)-\beta}{\Lambda}n+\sum_i\lambda_ic_i,\qquad
\frac{dc_i}{dt}=\frac{\beta_i}{\Lambda}n-\lambda_ic_i$$

6 gecikmiş grup (Keepin U-235 seti, β=650 pcm), Λ kullanıcı girdisi
(PWR ~2×10⁻⁵ s). **Sayısal yöntem:** çıktı adımı boyunca ρ donuk kabul
edilip 7×7 sistemin **kesin çözümü** y←e^{AΔt}y özayrışımla uygulanır —
prompt modun sertliği (stiffness) sorun olmaz. Doppler geri beslemesi:
ρ = ρ_dış + α_D(T−T₀); toplu yakıt ısısı
M·cp·dT/dt = (P−P₀) − M·cp(T−T₀)/τ_c.

**Doğrulama:** sabit ρ için benzetimin asimptotik periyodu analitik
**inhour** kökleriyle karşılaştırılır (verify.py: ±100 pcm'de sapma
<%1; −100 pcm için asimptota ulaşmak ~900 s gerektirir — test buna göre
kurgulanmıştır). Örnek ölçümler: +100 pcm → T=54.9 s (inhour 54.9);
+300 pcm fırlatma + Doppler(−2.5 pcm/K) @160 MW → tepe ≈ 316 MW,
kendiliğinden sınırlanma (+37 pcm net). Bozunum ısısı MODELLENMEZ;
scram sonrası güç yalnız kinetik seviyesidir.

## 5 · Mutlak birimler (işletme noktası)

Verilen kor termal gücünden: q‴_ort = P/V_yakıt; hücre bazında
q‴ = q‴_ort·(göreli güç). Akı ölçeği κ=200 MeV/fisyon ve ν̄=2.43
varsayımıyla: φ_abs = φ_rel·q‴_ort·ν̄/κ. Tipik PWR sağlaması:
q‴ ~ 50–110 W/cm³, termal akı ~3–6×10¹³ n/cm²·s.

## 6 · Sınırlar (dürüst kapsam beyanı)

1. **Difüzyon** teorisi: pin-ölçeği güçlü heterojenlikte yerel güçler
   sapar — C5G7 demosu bu sınırı ÖLÇER (k'da +31 pcm model hatası,
   pin gücünde ~%3; transport için ayrı S_N çözücümüz mevcuttur).
2. Donmuş 2-grup spektrum (tükenimde spektral kayma ihmal).
3. T-H geri beslemesi yok (transient'teki toplu Doppler modeli hariç);
   basınç/DNB hesabı kapsam dışıdır.
4. Bozunum ısısı ve uzaysal kinetik yok (nokta yaklaşımı).
5. Tükenim bu sürümde radyaldir (2-D); Gd/IFBA ayrık zehirler yok.

Bu sınırların her biri ya bir doğrulama vakasıyla nicelendirilmiştir ya
da açıkça kapsam dışı ilan edilmiştir; hiçbir sonuç sınırın ötesinde
bir doğruluk iddiasıyla sunulmaz.
