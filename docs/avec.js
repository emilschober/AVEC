// AVEC shared page script, loaded at the top of <body> on every page: applies the light/dark
// theme before the page paints (remembered in localStorage 'theme', otherwise following the
// system setting) and marks the current page in the navigation.
(function () {
    var KEY = 'theme';

    function stored() {
        try { return localStorage.getItem(KEY); } catch (e) { return null; }
    }

    function apply(theme) {
        var dark = theme === 'dark';
        document.body.classList.toggle('dark', dark);
        var btn = document.getElementById('theme-toggle');
        if (btn) btn.setAttribute('aria-pressed', String(dark));
    }

    var theme = stored();
    if (!theme) theme = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
    apply(theme);

    document.addEventListener('DOMContentLoaded', function () {
        apply(document.body.classList.contains('dark') ? 'dark' : 'light');
        var btn = document.getElementById('theme-toggle');
        if (btn) {
            btn.addEventListener('click', function () {
                var next = document.body.classList.contains('dark') ? 'light' : 'dark';
                try { localStorage.setItem(KEY, next); } catch (e) { /* theme just isn't remembered */ }
                var reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
                if (!document.startViewTransition || reduce) {
                    apply(next);
                    return;
                }
                // The new theme opens out from the toggle as a growing box
                var box = btn.getBoundingClientRect();
                var x = box.left + box.width / 2, y = box.top + box.height / 2;
                var start = 'inset(' + y + 'px ' + (window.innerWidth - x) + 'px ' + (window.innerHeight - y) + 'px ' + x + 'px)';
                var root = document.documentElement;
                root.classList.add('theme-switching');
                var transition = document.startViewTransition(function () { apply(next); });
                var done = function () { root.classList.remove('theme-switching'); };
                // A transition the browser skips (e.g. the tab is hidden) still applies the theme
                transition.ready.then(function () {
                    root.animate(
                        { clipPath: [start, 'inset(0px 0px 0px 0px)'] },
                        { duration: 550, easing: 'cubic-bezier(.2,.7,.2,1)', pseudoElement: '::view-transition-new(root)' });
                }).catch(done);
                transition.finished.then(done, done);
            });
        }
        var path = location.pathname.replace(/\/+$/, '') || '/';
        document.querySelectorAll('.site-nav a').forEach(function (a) {
            if (a.getAttribute('href') === path) a.setAttribute('aria-current', 'page');
        });

        // Tab bars (About page): an underline that slides to the active tab
        var tabs = document.querySelector('.tab-nav');
        if (tabs) {
            var line = document.createElement('span');
            line.className = 'tab-indicator';
            line.setAttribute('aria-hidden', 'true');
            tabs.appendChild(line);
            tabs.classList.add('has-indicator');
            var place = function () {
                var tab = tabs.querySelector('.tab-btn.active');
                if (!tab) return;
                line.style.width = tab.offsetWidth + 'px';
                line.style.transform = 'translate(' + tab.offsetLeft + 'px,' + (tab.offsetTop + tab.offsetHeight - 2) + 'px)';
            };
            place();
            // Slides only after the first placement
            requestAnimationFrame(function () { line.classList.add('is-ready'); });
            new MutationObserver(place).observe(tabs, { subtree: true, attributes: true, attributeFilter: ['class'] });
            window.addEventListener('resize', place);
        }
    });
})();
